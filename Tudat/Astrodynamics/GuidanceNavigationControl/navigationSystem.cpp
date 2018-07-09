#include <iostream>

#include <Eigen/LU>

#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/BasicMathematics/nonLinearLeastSquaresEstimation.h"
#include "Tudat/SimulationSetup/PropagationSetup/createEnvironmentUpdater.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to remove the error in gyroscope measurement based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromInertialMeasurementUnitMeasurement( const Eigen::Vector3d& currentInertialMeasurementUnitMeasurement,
                                                                    const Eigen::Vector6d& inertialMeasurementUnitErrors )
{
    return ( Eigen::Matrix3d::Identity( ) -
             Eigen::Matrix3d( inertialMeasurementUnitErrors.segment( 3, 3 ).asDiagonal( ) ) ) * // binomial approximation
            ( currentInertialMeasurementUnitMeasurement - inertialMeasurementUnitErrors.segment( 0, 3 ) );
}

//! Function to create navigation filter and root-finder objects for onboard state estimation.
void NavigationSystem::createNavigationSystemObjects(
        const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardSystemModel,
        const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardMeasurementModel )
{
    // Create filter object
    navigationFilter_ = filters::createFilter< >( navigationFilterSettings_, onboardSystemModel, onboardMeasurementModel );

    // Set initial time
    currentTime_ = navigationFilter_->getInitialTime( );
    currentOrbitCounter_ = 0;

    // Retrieve navigation filter step size and estimated state
    navigationRefreshStepSize_ = navigationFilter_->getIntegrationStepSize( );
    Eigen::Vector16d initialEstimatedState = navigationFilter_->getCurrentStateEstimate( );

    // Set initial rotational state
    currentEstimatedRotationalState_.setZero( );
    currentEstimatedRotationalState_.segment( 0, 4 ) = initialEstimatedState.segment( 6, 4 ).normalized( );

    // Set initial translational state
    setCurrentEstimatedCartesianState( initialEstimatedState.segment( 0, 6 ) );
    // this function also automatically stores the full state estimates at the current time

    // Update body and acceleration maps
    updateOnboardModel( );

    // Create root-finder object for bisection of aerodynamic acceleration curve
    // The values inserted are the tolerance in independent value (i.e., about twice the difference between
    // two time steps) and the maximum number of interations (i.e., 25 iterations)
    areaBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< > >(
                2.0 * navigationRefreshStepSize_ / currentTime_, 25 );
}

//! Function to create the onboard environment updater.
void NavigationSystem::createOnboardEnvironmentUpdater( )
{
    // Set integrated type and body list
    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedTypeAndBodyList;
    integratedTypeAndBodyList[ propagators::translational_state ] = { std::make_pair( spacecraftName_, "" ) };
    integratedTypeAndBodyList[ propagators::rotational_state ] = { std::make_pair( spacecraftName_, "" ) };

    // Create environment settings
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
            propagators::createTranslationalEquationsOfMotionEnvironmentUpdaterSettings( onboardAccelerationModelMap_, onboardBodyMap_ );

    // Create environment updater
    onboardEnvironmentUpdater_ = boost::make_shared< propagators::EnvironmentUpdater< double, double > >(
                onboardBodyMap_, environmentModelsToUpdate, integratedTypeAndBodyList );
}

//! Function to remove and calibrate (first time only) accelerometer errors.
void NavigationSystem::postProcessAccelerometerMeasurements(
        std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface,
        const std::map< double, Eigen::Vector3d >& mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface,
        const std::map< double, Eigen::Vector7d >& mapOfEstimatedRotationalStatesBelowAtmosphericInterface )
{
    // Inform user
    std::cout << "Removing Accelerometer Errors." << std::endl;

    // Apply smoothing method to noisy accelerometer data
    for ( unsigned int i = 0; i < vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
    {
        std::cout << vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ).transpose( ) << std::endl;
    }

    // Calibrate accelerometer errors if it is the first orbit
    if ( !atmosphereEstimatorInitialized_ )
    {
        // Inform user
        std::cout << "Calibrating Accelerometer" << std::endl;

        // Initial estimate on error values
        Eigen::Vector6d initialErrorEstimate = Eigen::Vector6d::Zero( );

        // Use non-linear least squares to solve for optimal value of errors
        estimatedAccelerometerErrors_ =
                linear_algebra::nonLinearLeastSquaresFit( boost::bind( &accelerometerErrorEstimationFunction, _1,
                                                                       vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface ),
                                                          initialErrorEstimate,
                                                          utilities::createConcatenatedEigenMatrixFromMapValues< double, double, 3 >(
                                                              mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface ) );
    }

    // Remove errors from accelerometer measurements and convert to inertial frame
    unsigned int i = 0;
    for ( rotationalStateIterator_ = mapOfEstimatedRotationalStatesBelowAtmosphericInterface.begin( );
          rotationalStateIterator_ != mapOfEstimatedRotationalStatesBelowAtmosphericInterface.end( ); rotationalStateIterator_++, i++ )
    {
        vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ) =
                removeErrorsFromInertialMeasurementUnitMeasurement( vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ),
                                                                    estimatedAccelerometerErrors_ );
        vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ) =
                linear_algebra::convertVectorToQuaternionFormat(
                    rotationalStateIterator_->second.segment( 0, 4 ) ).toRotationMatrix( ).transpose( ) *
                vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i );
    }
}

//! Function to run the Periapse Time Estimator (PTE).
void NavigationSystem::runPeriapseTimeEstimator(
        const std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
        const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface )
{
    // Inform user
    std::cout << "Running Periapse Time Estimator." << std::endl;

    // Separate time and accelerations
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > pairOfEstimatedKeplerianStateBelowAtmosphericInterface =
            utilities::extractKeyAndValuesFromMap< double, double, 6 >( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface );
    Eigen::VectorXd timesBelowAtmosphericInterface = pairOfEstimatedKeplerianStateBelowAtmosphericInterface.first;
    std::vector< double > vectorOfTimesBelowAtmosphericInterface =
            utilities::convertEigenVectorToStlVector( timesBelowAtmosphericInterface );
    Eigen::MatrixXd estimatedKeplerianStateBelowAtmosphericInterface = pairOfEstimatedKeplerianStateBelowAtmosphericInterface.second;
    Eigen::VectorXd estimatedTrueAnomalyBelowAtmosphericInterface = estimatedKeplerianStateBelowAtmosphericInterface.row( 5 ).transpose( );

    // Check that true anomaly and aerodynamic acceleration have the same length
    if ( vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.size( ) !=
         static_cast< unsigned int >( estimatedTrueAnomalyBelowAtmosphericInterface.rows( ) ) )
    {
        throw std::runtime_error( "Error in periapse time estimator. The sizes of the true anomaly and aerodynamic accelerations "
                                  "vectors do not match." );
    }

    // Set root-finder boundaries as the first and last times
    areaBisectionRootFinder_->resetBoundaries(
                vectorOfTimesBelowAtmosphericInterface.front( ), vectorOfTimesBelowAtmosphericInterface.back( ) );

    // Get intermediate variables
    Eigen::Vector6d initialEstimatedKeplerianState = estimatedKeplerianStateBelowAtmosphericInterface.col( 0 );

    // Set root-finder function as the area below the acceleration curve
    double estimatedActualPeriapseTime = areaBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< > >(
                    boost::bind( &areaBisectionFunction, _1, navigationRefreshStepSize_, timesBelowAtmosphericInterface,
                                 vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ) ) );

    // Interpolate to find estimated error in true anomaly
    double estimatedErrorInTrueAnomaly = interpolators::CubicSplineInterpolator< double, double >(
                vectorOfTimesBelowAtmosphericInterface,
                utilities::convertEigenVectorToStlVector( estimatedTrueAnomalyBelowAtmosphericInterface ) ).interpolate(
                estimatedActualPeriapseTime );
    std::cout << "Estimated Error in True Anomaly: " <<
                 unit_conversions::convertRadiansToDegrees( estimatedErrorInTrueAnomaly ) << " deg" << std::endl;
    // note that this represents directly the error in estimated true anomaly, since the true anomaly of
    // periapsis is zero by definition

    // Compute estimated change in velocity (i.e., Delta V) due to aerodynamic acceleration
    double estimatedChangeInVelocity = - numerical_quadrature::performExtendedSimpsonsQuadrature(
                navigationRefreshStepSize_, vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
    std::cout << "Estimated Change in Velocity: " << estimatedChangeInVelocity << " m/s" << std::endl;

    // Compute estimated mean motion by using the semi-major axis at beginning of atmospheric phase
    double initialEstimatedMeanMotion = std::sqrt( planetaryGravitationalParameter_ /
                                                   std::pow( initialEstimatedKeplerianState[ 0 ], 3 ) );

    // Compute estimated change in semi-major axis due to change in velocity
    // This value is computed by assuming that the change in velocity occurs at pericenter and is given by
    // by an impulsive shot. Also here, the eccentricity at the beginning of the atmospheric phase is used
    double estimatedChangeInSemiMajorAxisDueToChangeInVelocity = 2.0 / initialEstimatedMeanMotion * std::sqrt(
                ( 1.0 + initialEstimatedKeplerianState[ 1 ] ) /
            ( 1.0 - initialEstimatedKeplerianState[ 1 ] ) ) * estimatedChangeInVelocity;
    std::cout << "Estimated Change in Semi-major Axis: " <<
                 estimatedChangeInSemiMajorAxisDueToChangeInVelocity << " m" << std::endl;

    // Get estimated chagne in semi-major axis from estimated Keplerian state and find estimated error in semi-major axis
    double estimatedChangeInSemiMajorAxisFromKeplerianStateHistory = initialEstimatedKeplerianState[ 0 ] -
            estimatedKeplerianStateBelowAtmosphericInterface( 0, estimatedKeplerianStateBelowAtmosphericInterface.cols( ) - 1 );
    double estimatedErrorInSemiMajorAxis = estimatedChangeInSemiMajorAxisFromKeplerianStateHistory -
            estimatedChangeInSemiMajorAxisDueToChangeInVelocity;
    std::cout << "Estimated Error in Semi-major Axis: " << estimatedErrorInSemiMajorAxis << " m" << std::endl;

    // Compute estimated change in eccentricity due to change in velocity
    // The same assumption as for the case above holds
    double estimatedChangeInEccentricityDueToChangeInVelocity = 2.0 * std::sqrt( initialEstimatedKeplerianState[ 0 ] *
            ( 1 - std::pow( initialEstimatedKeplerianState[ 1 ], 2 ) ) / planetaryGravitationalParameter_ ) *
            estimatedChangeInVelocity;
    std::cout << "Estimated Change in Eccentricity: " <<
                 estimatedChangeInEccentricityDueToChangeInVelocity << std::endl;

    // Get estimated chagne in semi-major axis from estimated Keplerian state and find estimated error in semi-major axis
    double estimatedChangeInEccentricityFromKeplerianStateHistory = initialEstimatedKeplerianState[ 1 ] -
            estimatedKeplerianStateBelowAtmosphericInterface( 1, estimatedKeplerianStateBelowAtmosphericInterface.cols( ) - 1 );
    double estimatedErrorInEccentricity = estimatedChangeInEccentricityFromKeplerianStateHistory -
            estimatedChangeInEccentricityDueToChangeInVelocity;
    std::cout << "Estimated Error in Eccentricity: " << estimatedErrorInEccentricity << std::endl;

    // Combine errors to produce a vector of estimated error in Keplerian state
    Eigen::Vector6d estimatedErrorInKeplerianState = Eigen::Vector6d::Zero( );
    estimatedErrorInKeplerianState[ 0 ] = estimatedErrorInSemiMajorAxis;
    estimatedErrorInKeplerianState[ 1 ] = estimatedErrorInEccentricity;
    estimatedErrorInKeplerianState[ 5 ] = estimatedErrorInTrueAnomaly;
    historyOfEstimatedErrorsInKeplerianState_[ currentOrbitCounter_ ] = estimatedErrorInKeplerianState;

    // Compute conditions for propagation
    double timeFromPeriapsisToCurrentPosition = vectorOfTimesBelowAtmosphericInterface.back( ) - estimatedActualPeriapseTime;
    initialEstimatedKeplerianState[ 5 ] = 0.0;
    // set true anomaly to periapsis, to assume that no change in orbital elements has occurred during atmospheric phase
    Eigen::Vector6d correctedInitialEstimatedKeplerianState = initialEstimatedKeplerianState + estimatedErrorInKeplerianState;

    // Correct latest estimated Keplerian state with new information from PTE
    Eigen::Vector6d initialEstimatedCartesianElements =
            orbital_element_conversions::convertKeplerianToCartesianElements( initialEstimatedKeplerianState,
                                                                              planetaryGravitationalParameter_ );
    Eigen::Vector6d updatedCurrentEstimatedCartesianState =
            orbital_element_conversions::convertKeplerianToCartesianElements( correctedInitialEstimatedKeplerianState,
                                                                              planetaryGravitationalParameter_ );
    Eigen::Vector6d estimatedErrorInCartesianState = updatedCurrentEstimatedCartesianState - initialEstimatedCartesianElements;
    double estimatedAltitude = initialEstimatedCartesianElements.segment( 0, 3 ).norm( ) - planetaryRadius_;
    std::cout << updatedCurrentEstimatedCartesianState.transpose( ) << std::endl;
    std::cout << onboardBodyMap_.at( planetName_ )->getAtmosphereModel( )->getDensity( estimatedAltitude, 0.0, 0.0, 0.0 ) << std::endl;
    updatedCurrentEstimatedCartesianState += stateTransitionMatrixFunction_(
                updatedCurrentEstimatedCartesianState,
                onboardBodyMap_.at( planetName_ )->getAtmosphereModel( )->getDensity( estimatedAltitude, 0.0, 0.0, 0.0 ) ) *
            estimatedErrorInCartesianState * timeFromPeriapsisToCurrentPosition;

    // Update navigation system state estimates
    setCurrentEstimatedCartesianState( updatedCurrentEstimatedCartesianState );
    std::cout << "Before: " << initialEstimatedKeplerianState.transpose( ) << std::endl;
    std::cout << "Mid: " << correctedInitialEstimatedKeplerianState.transpose( ) << std::endl;
    std::cout << "After: " << currentEstimatedKeplerianState_.transpose( ) << std::endl;

    // Update navigation filter state and covariance
    Eigen::Vector16d updatedCurrentEstimatedState = navigationFilter_->getCurrentStateEstimate( );
    updatedCurrentEstimatedState.segment( 0, 6 ) = currentEstimatedCartesianState_;
    navigationFilter_->modifyCurrentStateAndCovarianceEstimates( updatedCurrentEstimatedState,
                                                                 Eigen::Matrix16d::Identity( ) );
    // the covariance matrix is reset to the identity, since the new state is improved in accuracy
}

//! Function to run the Atmosphere Estimator (AE).
void NavigationSystem::runAtmosphereEstimator(
        const std::map< double, Eigen::Vector6d >& mapOfEstimatedCartesianStatesBelowAtmosphericInterface,
        const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface )
{
    // Retrieve some physical parameters of the spacecraft
    double spacecraftMass = onboardBodyMap_.at( spacecraftName_ )->getBodyMass( );
    double referenceAerodynamicArea = onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
    double aerodynamicCoefficientsNorm =
            onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getCurrentForceCoefficients( ).norm( );

    // Pre-allocate variables
    std::vector< double > vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface;
    std::vector< double > vectorOfEstimatedAltitudesBelowAtmosphericInterface;

    // Convert estimated aerodynamic acceleration to estimated atmospheric density and compute altitude below atmospheric interface
    unsigned int i = 0;
    double currentRadialDistance;
    for ( std::map< double, Eigen::Vector6d >::const_iterator cartesianStateIterator =
          mapOfEstimatedCartesianStatesBelowAtmosphericInterface.begin( );
          cartesianStateIterator != mapOfEstimatedCartesianStatesBelowAtmosphericInterface.end( ); cartesianStateIterator++, i++ )
    {
        currentRadialDistance = cartesianStateIterator->second.segment( 0, 3 ).norm( );
        if ( currentRadialDistance < reducedAtmosphericInterfaceRadius_ )
        {
            // Get estimated density
            vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface.push_back(
                        2.0 * spacecraftMass / referenceAerodynamicArea / aerodynamicCoefficientsNorm /
                        cartesianStateIterator->second.segment( 3, 3 ).squaredNorm( ) *
                        vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.at( i ) );

            // Get estimated altitude
            vectorOfEstimatedAltitudesBelowAtmosphericInterface.push_back( currentRadialDistance - planetaryRadius_ );
        }
    }

    // Only proceed if satellite flew below reduced atmospheric interface altitude
    if ( vectorOfEstimatedAltitudesBelowAtmosphericInterface.size( ) > 0 )
    {
        // Convert vectors to Eigen
        Eigen::VectorXd estimatedAtmosphericDensitiesBelowAtmosphericInterface =
                utilities::convertStlVectorToEigenVector( vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface );
        Eigen::VectorXd estimatedAltitudesBelowAtmosphericInterface =
                utilities::convertStlVectorToEigenVector( vectorOfEstimatedAltitudesBelowAtmosphericInterface );

        // Find periapsis altitude
        double estimatedPeriapsisAltitude = estimatedAltitudesBelowAtmosphericInterface.minCoeff( );
        //        std::cout << "Periapsis: " << estimatedPeriapsisAltitude << std::endl;
        //        std::cout << "Altitudes: " << estimatedAltitudesBelowAtmosphericInterface.transpose( ) << std::endl;
        //        std::cout << "Densities: " << estimatedAtmosphericDensitiesBelowAtmosphericInterface.transpose( ) << std::endl;
        //        std::cout << "Densities: " << estimatedAtmosphericDensitiesBelowAtmosphericInterface.array( ).log( ) << std::endl;

        // Run least squares estimation process based on selected atmosphere model
        Eigen::VectorXd modelSpecificParameters;
        std::vector< double > vectorOfModelSpecificParameters;
        switch ( selectedOnboardAtmosphereModel_ )
        {
        case aerodynamics::exponential_atmosphere_model:
        case aerodynamics::three_wave_atmosphere_model:
        {
            // Compute information matrix
            Eigen::MatrixXd informationMatrix;
            informationMatrix.resize( estimatedAltitudesBelowAtmosphericInterface.rows( ), 2 );
            for ( unsigned int i = 0; i < estimatedAltitudesBelowAtmosphericInterface.rows( ); i++ )
            {
                informationMatrix.row( i ) << 1.0, estimatedPeriapsisAltitude - estimatedAltitudesBelowAtmosphericInterface[ i ];
            }

            // Use least squares polynomial fit
            Eigen::VectorXd logarithmOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface =
                    estimatedAtmosphericDensitiesBelowAtmosphericInterface.array( ).log( );
            Eigen::Vector2d estimatedAtmosphereModelParameters = ( informationMatrix.transpose( ) * informationMatrix ).inverse( ) *
                    informationMatrix.transpose( ) * logarithmOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface;

            // Add reference altitude to list of parameters and revert from logarithmic space
            if ( selectedOnboardAtmosphereModel_ == aerodynamics::exponential_atmosphere_model )
            {
                // If the exponential model is selected, add only the two estimated parameters
                modelSpecificParameters.resize( 3 );
                modelSpecificParameters[ 0 ] = estimatedPeriapsisAltitude;
                modelSpecificParameters[ 1 ] = std::exp( estimatedAtmosphereModelParameters[ 0 ] );
                modelSpecificParameters[ 2 ] = 1.0 / estimatedAtmosphereModelParameters[ 1 ];
            }
            else if ( selectedOnboardAtmosphereModel_ == aerodynamics::three_wave_atmosphere_model )
            {
                // If the three-wave model is selected, also add the extra two parameters
                modelSpecificParameters.resize( 5 );
                modelSpecificParameters[ 0 ] = estimatedPeriapsisAltitude;
                modelSpecificParameters[ 1 ] = std::exp( estimatedAtmosphereModelParameters[ 0 ] );
                modelSpecificParameters[ 2 ] = 1.0 / estimatedAtmosphereModelParameters[ 1 ];
                modelSpecificParameters[ 3 ] = 1.0;
                modelSpecificParameters[ 4 ] = 0.0;
            }
            break;
        }
        case aerodynamics::three_term_atmosphere_model:
        {
            // Initial estimate on atmosphere model parameters
            Eigen::Vector5d initialParameterEstimates;
            initialParameterEstimates << std::log( 2.424e-08 ), 1.0 / 6533.0, -1.0, 0.0, 0.0;

            // Use non-linear least squares to solve for optimal value of errors
            Eigen::Vector5d estimatedAtmosphereModelParameters = linear_algebra::nonLinearLeastSquaresFit(
                        boost::bind( &threeModelParametersEstimationFunction, _1, estimatedAltitudesBelowAtmosphericInterface,
                                     estimatedPeriapsisAltitude ),
                        initialParameterEstimates, estimatedAtmosphericDensitiesBelowAtmosphericInterface.array( ).log( ) );

            // Add reference altitude to list of parameters and revert from logarithmic space
            modelSpecificParameters.resize( 6 );
            modelSpecificParameters[ 0 ] = estimatedPeriapsisAltitude;
            modelSpecificParameters[ 1 ] = std::exp( estimatedAtmosphereModelParameters[ 0 ] );
            modelSpecificParameters[ 2 ] = 1.0 / estimatedAtmosphereModelParameters[ 1 ];
            modelSpecificParameters.segment( 3, 3 ) = estimatedAtmosphereModelParameters.segment( 2, 3 );
            break;
        }
        default:
            throw std::runtime_error( "Error in atmoshere estimation of navigation system. Atmosphere model not recognized." );
        }
        std::cout << "Atmosphere values: " << modelSpecificParameters.transpose( ) << std::endl;

        // Add values to hisotry
        historyOfEstimatedAtmosphereParameters_.push_back( modelSpecificParameters );

        // Perform moving average if enough parameters are available
        unsigned int numberOfSamplesForMovingAverage = historyOfEstimatedAtmosphereParameters_.size( );
        if ( numberOfSamplesForMovingAverage > 1 )
        {
            // If the number is larger than the number of samples to be used, use limiting value
            if ( numberOfSamplesForMovingAverage > numberOfAtmosphereSamplesForEstimation_ )
            {
                atmosphereEstimatorInitialized_ = true;
                numberOfSamplesForMovingAverage = numberOfAtmosphereSamplesForEstimation_;
            }

            // Compute moving average
            for ( unsigned int i = 1; i < numberOfSamplesForMovingAverage; i++ )
            {
                modelSpecificParameters += historyOfEstimatedAtmosphereParameters_.at(
                            historyOfEstimatedAtmosphereParameters_.size( ) - ( i + 1 ) );
            }
            modelSpecificParameters /= numberOfAtmosphereSamplesForEstimation_;
        }
        vectorOfModelSpecificParameters = utilities::convertEigenVectorToStlVector( modelSpecificParameters );

        // Update atmosphere settings of onboard body map
        onboardBodyMap_.at( planetName_ )->setAtmosphereModel(
                    boost::make_shared< aerodynamics::CustomConstantTemperatureAtmosphere >(
                        selectedOnboardAtmosphereModel_, 215.0, 197.0, 1.3, vectorOfModelSpecificParameters ) );
    }
}

} // namespace guidance_navigation_control

} // namespace tudat
