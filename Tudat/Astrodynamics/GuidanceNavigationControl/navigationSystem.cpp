#include <iostream>

#include <Eigen/LU>

#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"
#include "Tudat/Mathematics/BasicMathematics/nonLinearLeastSquaresEstimation.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/SimulationSetup/PropagationSetup/createEnvironmentUpdater.h"

#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to remove errors in inertial measurement unit measurements based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromInertialMeasurementUnitMeasurement( const Eigen::Vector3d& currentInertialMeasurementUnitMeasurement,
                                                                    const Eigen::Vector6d& inertialMeasurementUnitErrors )
{
    //    return ( Eigen::Matrix3d::Identity( ) -
    //             Eigen::Matrix3d( inertialMeasurementUnitErrors.segment( 3, 3 ).asDiagonal( ) ) ) * // binomial approximation
    //            ( currentInertialMeasurementUnitMeasurement - inertialMeasurementUnitErrors.segment( 0, 3 ) );
    return ( currentInertialMeasurementUnitMeasurement - inertialMeasurementUnitErrors.segment( 0, 3 ) );
}

//! Function to remove errors in the altimeter measurements based on the current orientation.
void removeErrorsFromAltimeterMeasurement( double& currentMeasuredAltitude, const double planetaryRadius )
{
    // Assume perfect pointing
    double currentPointingAngle = 0.0;

    // If pointing angle is larger than the singularity
    if ( std::fabs( currentPointingAngle ) > std::numeric_limits< double >::epsilon( ) )
    {
        // Remove pointing error from altimeter
        double sineOfCurrentPointingAngle = std::sin( currentPointingAngle );
        currentMeasuredAltitude = planetaryRadius * (
                    std::sin( mathematical_constants::PI - currentPointingAngle -
                              std::asin( currentMeasuredAltitude / planetaryRadius *
                                         sineOfCurrentPointingAngle ) ) / sineOfCurrentPointingAngle - 1.0 );
    }
}

//! Function to create navigation objects for onboard state estimation.
void NavigationSystem::createNavigationSystemObjects(
        const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardSystemModel,
        const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardMeasurementModel )
{
    // Create filter object
    switch ( navigationFilterSettings_->filteringTechnique_ )
    {
    case filters::extended_kalman_filter:
    {
        // Create inputs
        boost::function< double( const Eigen::Vector6d& ) > densityFunction =
                boost::bind( &NavigationSystem::getDensityAtSpecifiedConditions, this, _1 );
        double aerodynamicParameters = 0.5 *
                onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getCurrentForceCoefficients( )[ 0 ] *
                onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( ) /
                onboardBodyMap_.at( spacecraftName_ )->getBodyMass( );

        // Create filter
        navigationFilter_ = filters::createFilter< double, double >(
                    navigationFilterSettings_, onboardSystemModel, onboardMeasurementModel,
                    boost::bind( &computeSystemJacobianMatrix, _1, _2, densityFunction, planetaryGravitationalParameter_,
                                 planetaryRadius_, secondDegreeGravitationalMoment_, aerodynamicParameters ),
                    boost::lambda::constant( Eigen::Matrix12d::Identity( ) ),
                    boost::bind( &computeMeasurementJacobianMatrix, _1, _2, densityFunction, aerodynamicParameters ),
                    boost::lambda::constant( Eigen::Matrix3d::Identity( ) ) );
        break;
    }
    case filters::unscented_kalman_filter:
    {
        navigationFilter_ = filters::createFilter< double, double >( navigationFilterSettings_, onboardSystemModel, onboardMeasurementModel );
        break;
    }
    default:
        throw std::runtime_error( "Error in setting up navigation system. The requested filtering technique is not supported." );
    }

    // Set initial time
    currentTime_ = navigationFilter_->getInitialTime( );
    currentOrbitCounter_ = 0;

    // Retrieve navigation filter step size and estimated state
    navigationRefreshStepSize_ = navigationFilter_->getIntegrationStepSize( );
    Eigen::Vector12d initialEstimatedState = navigationFilter_->getCurrentStateEstimate( );

    // Set initial translational state
    setCurrentEstimatedCartesianState( initialEstimatedState.segment( 0, 6 ) );
    // this function also automatically stores the full state estimates at the current time

    // Store the apoapsis value of Keplerian state for the Periapse Time Estimator
    estimatedKeplerianStateAtPreviousApoapsis_ = currentEstimatedKeplerianState_;

    // Update body and acceleration maps
    updateOnboardModel( );

    // Create root-finder object for bisection of aerodynamic acceleration curve
    // The values inserted are the tolerance in independent value (i.e., about twice the difference between
    // two time steps) and the maximum number of iterations (i.e., 25 iterations)
    areaBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >(
                2.0 * navigationRefreshStepSize_ / currentTime_, 25 );

    // Create object of dependent variables to save
    std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                          propagators::local_dynamic_pressure_dependent_variable, spacecraftName_, planetName_ ) );
    dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                          propagators::local_aerodynamic_heat_rate_dependent_variable, spacecraftName_, planetName_ ) );

    // Create object for propagation of spacecraft state with user-provided initial conditions
    onboardIntegratorSettings_ =
            boost::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettings< double, Eigen::VectorXd > >(
                0.0, 10.0, numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg56, 1e-5, 1e5, 1.0e-12, 1.0e-12 );
    onboardPropagatorSettings_ = boost::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                std::vector< std::string >( 1, planetName_ ), onboardAccelerationModelMap_,
                std::vector< std::string >( 1, spacecraftName_ ), Eigen::Vector6d::Zero( ),
                0.0, propagators::unified_state_model_exponential_map,
                boost::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList, false ) );
}

//! Function to create the onboard environment updater.
void NavigationSystem::createOnboardEnvironmentUpdater( )
{
    // Set integrated type and body list
    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedTypeAndBodyList;
    integratedTypeAndBodyList[ propagators::translational_state ] = { std::make_pair( spacecraftName_, "" ) };

    // Create environment settings
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
            propagators::createTranslationalEquationsOfMotionEnvironmentUpdaterSettings( onboardAccelerationModelMap_, onboardBodyMap_ );

    // Create environment updater
    onboardEnvironmentUpdater_ = boost::make_shared< propagators::EnvironmentUpdater< double, double > >(
                onboardBodyMap_, environmentModelsToUpdate, integratedTypeAndBodyList );
}

//! Function to post-process the accelerometer measurements.
void NavigationSystem::postProcessAccelerometerMeasurements(
        std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface )
{
    // Inform user
    std::cout << std::endl << "Removing Accelerometer Errors." << std::endl;

    // Remove errors from accelerometer measurements and convert to inertial frame
    for ( unsigned int i = 0; i < vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
    {
        vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ) =
                removeErrorsFromInertialMeasurementUnitMeasurement(
                    vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ),
                    navigationFilter_->getCurrentStateEstimate( ).segment( 6, 6 ) );
    }

    // Apply smoothing method to noisy accelerometer data
    unsigned int numberOfSamplePoints = static_cast< unsigned int >( 50.0 / navigationRefreshStepSize_ );
    numberOfSamplePoints = ( ( numberOfSamplePoints % 2 ) == 0 ) ? numberOfSamplePoints + 1 : numberOfSamplePoints; // only odd values
    vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface =
            statistics::computeMovingAverage( vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface, numberOfSamplePoints );
}

//! Function to run the Periapse Time Estimator (PTE).
void NavigationSystem::runPeriapseTimeEstimator(
        std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
        const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface )
{
    // Inform user
    std::cout << std::endl << "Running Periapse Time Estimator." << std::endl;

    input_output::writeDataMapToTextFile( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                          "kepler_" + std::to_string( currentOrbitCounter_ ) + ".dat", "/Users/Michele/Desktop/" );
    input_output::writeMatrixToFile(
                utilities::convertStlVectorToEigenVector( vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ),
                "aero_" + std::to_string( currentOrbitCounter_ ) + ".dat", 16, "/Users/Michele/Desktop/" );

    // Separate time and accelerations
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > pairOfEstimatedKeplerianStateBelowAtmosphericInterface =
            utilities::extractKeyAndValuesFromMap< double, double, 6 >( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface );
    Eigen::VectorXd timesBelowAtmosphericInterface = pairOfEstimatedKeplerianStateBelowAtmosphericInterface.first;
    std::vector< double > vectorOfTimesBelowAtmosphericInterface = utilities::convertEigenVectorToStlVector( timesBelowAtmosphericInterface );
    Eigen::MatrixXd estimatedKeplerianStateBelowAtmosphericInterface = pairOfEstimatedKeplerianStateBelowAtmosphericInterface.second;
    Eigen::VectorXd estimatedTrueAnomalyBelowAtmosphericInterface = estimatedKeplerianStateBelowAtmosphericInterface.row( 5 ).transpose( );

    // Check that true anomaly and aerodynamic acceleration have the same length
    if ( vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.size( ) !=
         static_cast< unsigned int >( estimatedTrueAnomalyBelowAtmosphericInterface.rows( ) ) )
    {
        throw std::runtime_error( "Error in Periapse Time Estimator. The sizes of the true anomaly and aerodynamic acceleration "
                                  "vectors do not match." );
    }

    // Set root-finder boundaries as the first and last times
    areaBisectionRootFinder_->resetBoundaries(
                vectorOfTimesBelowAtmosphericInterface.front( ), vectorOfTimesBelowAtmosphericInterface.back( ) );

    // Set root-finder function as the area below the acceleration curve
    double estimatedActualPeriapseTime = areaBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
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

    // Propagate state to apoapsis to see what will be the value of semi-major axis and eccentricity
    boost::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::keplerian_state_dependent_variable, spacecraftName_, planetName_, 5 ),
                mathematical_constants::PI, false );//, true,
//                boost::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 0.001, 25 ) );
    Eigen::Vector6d estimatedCartesianStateAtNextApoapsis =
            propagateTranslationalStateWithCustomTerminationSettings( terminationSettings ).first.rbegin( )->second;
    Eigen::Vector6d estimatedKeplerianStateAtNextApoapsis =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                estimatedCartesianStateAtNextApoapsis, planetaryGravitationalParameter_ );
    std::cout << "Propagated: " << estimatedKeplerianStateAtNextApoapsis.transpose( ) << std::endl
              << "Initial: " << currentEstimatedKeplerianState_.transpose( ) << std::endl
              << "Previous: " << estimatedKeplerianStateAtPreviousApoapsis_.transpose( ) << std::endl;

    // Compute estimated change in velocity (i.e., Delta V) due to aerodynamic acceleration
    double estimatedChangeInVelocity = - periapseEstimatorConstant_ * numerical_quadrature::performExtendedSimpsonsQuadrature(
                navigationRefreshStepSize_, vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
    std::cout << "Estimated Change in Velocity: " << estimatedChangeInVelocity << " m/s" << std::endl;

    // Compute estimated mean motion by using the semi-major axis at beginning of atmospheric phase
    double initialEstimatedMeanMotion = std::sqrt( planetaryGravitationalParameter_ /
                                                   std::pow( estimatedKeplerianStateAtPreviousApoapsis_[ 0 ], 3 ) );

    // Compute estimated change in semi-major axis due to change in velocity
    // This value is computed by assuming that the change in velocity occurs at pericenter and is given by
    // by an impulsive shot. Also here, the eccentricity at the beginning of the atmospheric phase is used
    double estimatedChangeInSemiMajorAxisDueToChangeInVelocity = 2.0 / initialEstimatedMeanMotion * std::sqrt(
                ( 1.0 + estimatedKeplerianStateAtPreviousApoapsis_[ 1 ] ) /
            ( 1.0 - estimatedKeplerianStateAtPreviousApoapsis_[ 1 ] ) ) * estimatedChangeInVelocity;
    std::cout << "Estimated Change in Semi-major Axis: " <<
                 estimatedChangeInSemiMajorAxisDueToChangeInVelocity << " m" << std::endl;

    // Get estimated chagne in semi-major axis from estimated Keplerian state and find estimated error in semi-major axis
    double estimatedChangeInSemiMajorAxisFromKeplerianStateHistory =
            estimatedKeplerianStateAtNextApoapsis[ 0 ] - estimatedKeplerianStateAtPreviousApoapsis_[ 0 ];
    double estimatedErrorInSemiMajorAxis = estimatedChangeInSemiMajorAxisFromKeplerianStateHistory -
            estimatedChangeInSemiMajorAxisDueToChangeInVelocity;
    std::cout << "Estimated Error in Semi-major Axis: " << estimatedErrorInSemiMajorAxis << " m" << std::endl;

    // Compute estimated change in eccentricity due to change in velocity
    // The same assumption as for the case above holds
    double estimatedChangeInEccentricityDueToChangeInVelocity = 2.0 * std::sqrt( estimatedKeplerianStateAtPreviousApoapsis_[ 0 ] *
            ( 1 - std::pow( estimatedKeplerianStateAtPreviousApoapsis_[ 1 ], 2 ) ) / planetaryGravitationalParameter_ ) *
            estimatedChangeInVelocity;
    std::cout << "Estimated Change in Eccentricity: " <<
                 estimatedChangeInEccentricityDueToChangeInVelocity << std::endl;

    // Get estimated chagne in semi-major axis from estimated Keplerian state and find estimated error in semi-major axis
    double estimatedChangeInEccentricityFromKeplerianStateHistory =
            estimatedKeplerianStateAtNextApoapsis[ 1 ] - estimatedKeplerianStateAtPreviousApoapsis_[ 1 ];
    double estimatedErrorInEccentricity = estimatedChangeInEccentricityFromKeplerianStateHistory -
            estimatedChangeInEccentricityDueToChangeInVelocity;
    std::cout << "Estimated Error in Eccentricity: " << estimatedErrorInEccentricity << std::endl;

    // Combine errors to produce a vector of estimated error in Keplerian state
    Eigen::Vector6d estimatedErrorInKeplerianState = Eigen::Vector6d::Zero( );
    estimatedErrorInKeplerianState[ 0 ] = estimatedErrorInSemiMajorAxis;//estimatedChangeInSemiMajorAxisDueToChangeInVelocity
    estimatedErrorInKeplerianState[ 1 ] = estimatedErrorInEccentricity;//estimatedChangeInEccentricityDueToChangeInVelocity
//    estimatedErrorInKeplerianState[ 5 ] = estimatedErrorInTrueAnomaly;
    historyOfEstimatedErrorsInKeplerianState_[ currentOrbitCounter_ ] = estimatedErrorInKeplerianState;

    // Compute updated estimate in Keplerian state at current time by removing the estimated error
    Eigen::Vector6d updatedCurrentKeplerianState = currentEstimatedKeplerianState_;
    updatedCurrentKeplerianState -= estimatedErrorInKeplerianState;
    // the updated state is initially defined as the one with semi-major axis and eccentricity from the time before the
    // atmospheric pass and inclination, right ascension of ascending node, argument of periapsis and true anomaly from the
    // latest estimate; then the error in semi-major axis, eccentricity and true anomaly is subtracted

    // Update navigation system state estimates
    Eigen::Matrix12d currentEstimatedCovarianceMatrix = navigationFilter_->getCurrentCovarianceEstimate( );
    currentEstimatedCovarianceMatrix.block( 0, 0, 6, 6 ).setIdentity( );
    setCurrentEstimatedKeplerianState( updatedCurrentKeplerianState, currentEstimatedCovarianceMatrix );
    // the covariance matrix is reset to the identity, since the new state is improved in accuracy

    // Correct history of Keplerian elements by removing error in true anomaly
    for ( std::map< double, Eigen::Vector6d >::iterator
          keplerianStateIterator = mapOfEstimatedKeplerianStatesBelowAtmosphericInterface.begin( );
          keplerianStateIterator != mapOfEstimatedKeplerianStatesBelowAtmosphericInterface.end( ); keplerianStateIterator++ )
    {
        keplerianStateIterator->second[ 5 ] -= estimatedErrorInTrueAnomaly;
    }
}

//! Function to run the Atmosphere Estimator (AE).
void NavigationSystem::runAtmosphereEstimator(
        const std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
        const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface )
{
    // Inform user
    std::cout << std::endl << "Estimating Atmospheric Parameters." << std::endl;

    // Retrieve some physical parameters of the spacecraft
    double spacecraftMass = onboardBodyMap_.at( spacecraftName_ )->getBodyMass( );
    double referenceAerodynamicArea = onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
    double dragCoefficient = onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getCurrentForceCoefficients( )[ 0 ];

    // Pre-allocate variables
    std::vector< double > vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface;
    std::vector< double > vectorOfEstimatedAltitudesBelowAtmosphericInterface;

    // Convert estimated aerodynamic acceleration to estimated atmospheric density and compute altitude below atmospheric interface
    unsigned int i = 0;
    double currentRadialDistance;
    for ( std::map< double, Eigen::Vector6d >::const_iterator
          keplerianStateIterator = mapOfEstimatedKeplerianStatesBelowAtmosphericInterface.begin( );
          keplerianStateIterator != mapOfEstimatedKeplerianStatesBelowAtmosphericInterface.end( ); keplerianStateIterator++, i++ )
    {
        currentRadialDistance = basic_astrodynamics::computeKeplerRadialDistance( keplerianStateIterator->second );
        if ( currentRadialDistance < reducedAtmosphericInterfaceRadius_ )
        {
            // Get estimated density
            vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface.push_back(
                        2.0 * spacecraftMass / referenceAerodynamicArea / dragCoefficient /
                        std::pow( basic_astrodynamics::computeKeplerOrbitalVelocity( keplerianStateIterator->second,
                                                                                     planetaryGravitationalParameter_ ), 2 ) *
                        vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.at( i ) );

            // Get estimated altitude
            vectorOfEstimatedAltitudesBelowAtmosphericInterface.push_back( currentRadialDistance - planetaryRadius_ );
        }
    }

    // Only proceed if satellite flew below reduced atmospheric interface altitude
    if ( !vectorOfEstimatedAltitudesBelowAtmosphericInterface.empty( ) )
    {
        // Convert vectors to Eigen
        Eigen::VectorXd estimatedAtmosphericDensitiesBelowAtmosphericInterface =
                utilities::convertStlVectorToEigenVector( vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface );
        Eigen::VectorXd estimatedAltitudesBelowAtmosphericInterface =
                utilities::convertStlVectorToEigenVector( vectorOfEstimatedAltitudesBelowAtmosphericInterface );

        // Find periapsis altitude
        double estimatedPeriapsisAltitude = estimatedAltitudesBelowAtmosphericInterface.minCoeff( );
//        std::cout << "Periapsis: " << estimatedPeriapsisAltitude << std::endl;
//        std::cout << "Altitudes: " << std::endl << estimatedAltitudesBelowAtmosphericInterface.transpose( ) << std::endl;
//        std::cout << "Altitudes: " << std::endl << estimatedAltitudesBelowAtmosphericInterface.array( ).pow( -2 ).transpose( ) << std::endl;
//        std::cout << "Densities: " << std::endl << estimatedAtmosphericDensitiesBelowAtmosphericInterface.transpose( ) << std::endl;
//        std::cout << "Densities: " << estimatedAtmosphericDensitiesBelowAtmosphericInterface.array( ).log( ) << std::endl;
//        std::cout << "Accelerations: " << std::endl << utilities::convertStlVectorToEigenVector(
//                         vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ).transpose( ) << std::endl;

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
            Eigen::Vector2d estimatedAtmosphereModelParameters = linear_algebra::performLeastSquaresAdjustmentFromInformationMatrix(
                        informationMatrix, logarithmOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface,
                        estimatedAltitudesBelowAtmosphericInterface.array( ).pow( -2 ), false ).first;

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
                modelSpecificParameters[ 3 ] = 1.0; // should be randomized
                modelSpecificParameters[ 4 ] = 0.0; // should depend on dust storms
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

        // Check that estimated atmospheric parameters make sense
        if ( ( modelSpecificParameters[ 1 ] < 0 ) || // density at zero altitude is negative
             ( modelSpecificParameters[ 2 ] < 0 ) ) // scale height is negative
        {
            // Inform user
            std::cerr << "Warning in atmosphere estimator. The estimated atmospheric parameters are not in line "
                         "with what is expected. Taking an average of the previous estimated values, instead." << std::endl;

            // Take average of all previous values
            modelSpecificParameters.setZero( );
            for ( unsigned int i = 0; i < historyOfEstimatedAtmosphereParameters_.size( ); i++ )
            {
                modelSpecificParameters += historyOfEstimatedAtmosphereParameters_.at( i );
            }
            modelSpecificParameters /= historyOfEstimatedAtmosphereParameters_.size( );
        }
        std::cout << "Atmosphere values: " << modelSpecificParameters.transpose( ) << std::endl;

        // Add values to hisotry
        historyOfEstimatedAtmosphereParameters_.push_back( modelSpecificParameters );

        // Perform moving average if enough parameters are available
        unsigned int numberOfSamplesForMovingAverage = historyOfEstimatedAtmosphereParameters_.size( );
        if ( numberOfSamplesForMovingAverage > 1 )
        {
            // If the number is larger than the number of samples to be used, use limiting value
            if ( numberOfSamplesForMovingAverage > numberOfRequiredAtmosphereSamplesForInitiation_ )
            {
                atmosphereEstimatorInitialized_ = true;
                numberOfSamplesForMovingAverage = numberOfRequiredAtmosphereSamplesForInitiation_;
            }

            // Compute moving average
            for ( unsigned int i = 1; i < numberOfSamplesForMovingAverage; i++ )
            {
                modelSpecificParameters += historyOfEstimatedAtmosphereParameters_.at(
                            historyOfEstimatedAtmosphereParameters_.size( ) - ( i + 1 ) );
            }
            modelSpecificParameters /= numberOfSamplesForMovingAverage;
        }
        vectorOfModelSpecificParameters = utilities::convertEigenVectorToStlVector( modelSpecificParameters );
        std::cout << "Averaged atmosphere values: " << modelSpecificParameters.transpose( ) << std::endl;

        // Update atmosphere settings of onboard body map
        onboardBodyMap_.at( planetName_ )->setAtmosphereModel(
                    boost::make_shared< aerodynamics::CustomConstantTemperatureAtmosphere >(
                        selectedOnboardAtmosphereModel_, 215.0, 197.0, 1.3, vectorOfModelSpecificParameters ) );
    }
}

} // namespace guidance_navigation_control

} // namespace tudat
