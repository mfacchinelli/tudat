#include <iostream>

#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"
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

//! Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
double areaBisectionFunction( const double currentTimeGuess, const double constantTimeStep,
                              const Eigen::VectorXd& onboardTime,
                              const std::vector< double >& estimatedAerodynamicAccelerationMagnitude )
{
    // Find nearest lower index to true anomaly guess
    int nearestLowerIndex = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch( onboardTime, currentTimeGuess );
    nearestLowerIndex = ( nearestLowerIndex == 0 ) ? 1 : nearestLowerIndex;

    // Compute trapezoidal quadrature to integrate the area until and after the current guess
    std::cout << "Guess: " << currentTimeGuess - 236304000.0 << std::endl
              << "Index: " << nearestLowerIndex << std::endl
              << "Total length: " << onboardTime.rows( ) << std::endl;
    double lowerSliceQuadratureResult = numerical_quadrature::performExtendedSimpsonsQuadrature(
                constantTimeStep, utilities::sliceStlVector( estimatedAerodynamicAccelerationMagnitude, 0, nearestLowerIndex ) );
    double upperSliceQuadratureResult = numerical_quadrature::performExtendedSimpsonsQuadrature(
                constantTimeStep, utilities::sliceStlVector( estimatedAerodynamicAccelerationMagnitude, nearestLowerIndex ) );
    std::cout << "Lower: " << lowerSliceQuadratureResult << std::endl;
    std::cout << "Upper: " << upperSliceQuadratureResult << std::endl;
    std::cout << "Result: " << upperSliceQuadratureResult - lowerSliceQuadratureResult << std::endl;

    // Return difference in areas
    return upperSliceQuadratureResult - lowerSliceQuadratureResult;
}

//! Function to be used as input to the non-linear least squares process to determine the accelerometer errors.
std::pair< Eigen::VectorXd, Eigen::MatrixXd > accelerometerErrorEstimationFunction(
        const Eigen::Vector6d& currentErrorEstimate,
        const std::vector< Eigen::Vector3d >& vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface )
{
    // Set variables to zero
    Eigen::VectorXd expectedAcceleration = Eigen::VectorXd::Zero(
                3 * vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.size( ) );
    Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Zero(
                3 * vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.size( ), 6 );

    // Loop over each acceleration to add values to matrix
    for ( unsigned int i = 0; i < vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
    {
        // Find current expected measurement
        expectedAcceleration.segment( 3 * i, 3 ) =
                ( Eigen::Matrix3d::Identity( ) - Eigen::Matrix3d( currentErrorEstimate.segment( 3, 3 ).asDiagonal( ) ) ) *
                ( vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.at( i ) - currentErrorEstimate.segment( 0, 3 ) );

        // Find current Jacobian matrix
        designMatrix.row( 3 * i ) << currentErrorEstimate[ 3 ] - 1.0, 0.0, 0.0, currentErrorEstimate[ 0 ] -
                vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.at( i )[ 0 ], 0.0, 0.0;
        designMatrix.row( 3 * i + 1 ) << 0.0, currentErrorEstimate[ 4 ] - 1.0, 0.0, 0.0, currentErrorEstimate[ 1 ] -
                vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.at( i )[ 1 ], 0.0;
        designMatrix.row( 3 * i + 2 ) << 0.0, 0.0, currentErrorEstimate[ 5 ] - 1.0, 0.0, 0.0, currentErrorEstimate[ 2 ] -
                vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.at( i )[ 2 ];
    }

    // Return acceleration and design matrix as pair
    return std::make_pair( expectedAcceleration, designMatrix );
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
    // two time steps) and the maximum number of interations (i.e., 10 iterations)
    areaBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >(
                2.0 * navigationRefreshStepSize_ / currentTime_, 10 );
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
void NavigationSystem::removeAccelerometerErrors(
        std::vector< Eigen::Vector3d >& vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface,
        const std::map< double, Eigen::Vector3d >& mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface )
{
    std::cout << "Removing Accelerometer Errors." << std::endl;

    // Calibrate accelerometer errors if it is the first orbit
    if ( currentOrbitCounter_ == 0 )
    {
        std::cout << "Calibrating Accelerometer" << std::endl;

        // Initial estimate on error values
        Eigen::Vector6d initialErrorEstimate = Eigen::Vector6d::Zero( );

        // Use non-linear least squares to solve for optimal value of errors
        accelerometerErrors_ =
                linear_algebra::nonLinearLeastSquaresFit( boost::bind( &accelerometerErrorEstimationFunction, _1,
                                                                       vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface ),
                                                          initialErrorEstimate,
                                                          utilities::createConcatenatedEigenMatrixFromMapValues< double, double, 3 >(
                                                              mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface ) );
    }

    // Remove errors from accelerometer measurements
    for ( std::vector< Eigen::Vector3d >::iterator
          vectorIterator = vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.begin( );
          vectorIterator != vectorOfEstimatedAerodynamicAccelerationBelowAtmosphericInterface.end( ); vectorIterator++ )
    {
        *vectorIterator = removeErrorsFromInertialMeasurementUnitMeasurement( *vectorIterator, accelerometerErrors_ );
    }
}

//! Function to run the Periapse Time Estimator (PTE).
void NavigationSystem::runPeriapseTimeEstimator(
        const std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
        const std::vector< double >& vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface )
{
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
    if ( vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.size( ) !=
         static_cast< unsigned int >( estimatedTrueAnomalyBelowAtmosphericInterface.rows( ) ) )
    {
        throw std::runtime_error( "Error in periapse time estimator. The sizes of the true anomaly and aerodynamic accelerations "
                                  "vectors do not match." );
    }

    // Set root-finder boundaries as the first and last times
    areaBisectionRootFinder_->resetBoundaries(
                vectorOfTimesBelowAtmosphericInterface.front( ), vectorOfTimesBelowAtmosphericInterface.back( ) );
    std::cout << "Init. time: " << vectorOfTimesBelowAtmosphericInterface.front( ) - 236304000.0 << std::endl
              << "Fin. time: " << vectorOfTimesBelowAtmosphericInterface.back( ) - 236304000.0 << std::endl;

    // Get intermediate variables
    Eigen::Vector6d initialEstimatedKeplerianState = estimatedKeplerianStateBelowAtmosphericInterface.col( 0 );

    // Set root-finder function as the area below the acceleration curve
    double estimatedActualPeriapseTime = areaBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< > >(
                    boost::bind( &areaBisectionFunction, _1, navigationRefreshStepSize_, timesBelowAtmosphericInterface,
                                 vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ) ) );

    // Interpolate to find estimated error in true anomaly
    double estimatedErrorInTrueAnomaly = interpolators::CubicSplineInterpolator< double, double >(
                vectorOfTimesBelowAtmosphericInterface,
                utilities::convertEigenVectorToStlVector( estimatedTrueAnomalyBelowAtmosphericInterface ) ).interpolate(
                estimatedActualPeriapseTime );
    std::cout << "Estimated Error in True Anomaly: " <<
                 unit_conversions::convertDegreesToRadians( estimatedErrorInTrueAnomaly ) << " deg" << std::endl;
    // note that this represents directly the error in estimated true anomaly, since the true anomaly of
    // periapsis is zero by definition

    // Find nearest lower index to error in true anomaly
    int estimatedPeriapsisIndex = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                estimatedTrueAnomalyBelowAtmosphericInterface, estimatedErrorInTrueAnomaly );
    std::cout << "Check: " << estimatedActualPeriapseTime - 236304000.0 << ", "
              << vectorOfTimesBelowAtmosphericInterface.at( estimatedPeriapsisIndex ) - 236304000.0 << std::endl;

    // Compute estimated change in velocity (i.e., Delta V) due to aerodynamic acceleration
    double estimatedChangeInVelocity = - numerical_quadrature::performExtendedSimpsonsQuadrature(
                navigationRefreshStepSize_, vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
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
            estimatedKeplerianStateBelowAtmosphericInterface( 0, estimatedKeplerianStateBelowAtmosphericInterface.cols( ) );
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
            estimatedKeplerianStateBelowAtmosphericInterface( 1, estimatedKeplerianStateBelowAtmosphericInterface.cols( ) );
    double estimatedErrorInEccentricity = estimatedChangeInEccentricityFromKeplerianStateHistory -
            estimatedChangeInEccentricityDueToChangeInVelocity;
    std::cout << "Estimated Error in Eccentricity: " << estimatedErrorInEccentricity << std::endl;

    // Combine errors to produce a vector of estimated error in Keplerian state
    Eigen::Vector6d estimatedErrorInKeplerianState = Eigen::Vector6d::Zero( );
    estimatedErrorInKeplerianState[ 0 ] = estimatedErrorInSemiMajorAxis;
    estimatedErrorInKeplerianState[ 1 ] = estimatedErrorInEccentricity;
    estimatedErrorInKeplerianState[ 5 ] = estimatedErrorInTrueAnomaly;

    // Correct latest estimated Keplerian state with new information from PTE
    double timeFromPeriapsisToCurrentPosition = vectorOfTimesBelowAtmosphericInterface.back( ) - estimatedActualPeriapseTime;
    Eigen::Vector6d updatedCurrentEstimatedKeplerianState = initialEstimatedKeplerianState;
    initialEstimatedKeplerianState[ 5 ] = 0.0; // set true anomaly to periapsis, to assume that no change in orbital elements
                                               // has occurred during atmospheric phase
    updatedCurrentEstimatedKeplerianState += stateTransitionMatrixFunction_( initialEstimatedKeplerianState ) *
            estimatedErrorInKeplerianState * timeFromPeriapsisToCurrentPosition;

    // Update navigation system state estimates
    setCurrentEstimatedKeplerianState( updatedCurrentEstimatedKeplerianState );
    historyOfEstimatedErrorsInKeplerianState_[ currentOrbitCounter_ ] = estimatedErrorInKeplerianState;

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
        const std::vector< double >& vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface )
{
    // Retrieve some physical parameters of the spacecraft
    double spacecraftMass = onboardBodyMap_.at( spacecraftName_ )->getBodyMass( );
    double referenceAerodynamicArea = onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getReferenceArea( );
    double aerodynamicCoefficientsNorm =
            onboardBodyMap_.at( spacecraftName_ )->getAerodynamicCoefficientInterface( )->getCurrentForceCoefficients( ).norm( );

    // Convert estimated aerodynamic acceleration to estimated atmospheric density and compute altitude below atmospheric interface
    unsigned int i = 0;
    std::vector< double > vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface;
    std::vector< double > vectorOfEstimatedAltitudesBelowAtmosphericInterface;
    for ( std::map< double, Eigen::Vector6d >::const_iterator cartesianStateIterator =
          mapOfEstimatedCartesianStatesBelowAtmosphericInterface.begin( );
          cartesianStateIterator != mapOfEstimatedCartesianStatesBelowAtmosphericInterface.end( ); cartesianStateIterator++, i++ )
    {
        // Get estimated density
        vectorOfEstimatedAtmosphericDensitiesBelowAtmosphericInterface.push_back(
                    2.0 * spacecraftMass / referenceAerodynamicArea / aerodynamicCoefficientsNorm /
                    cartesianStateIterator->second.segment( 3, 3 ).squaredNorm( ) *
                    vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.at( i ) );

        // Get estimated altitude
        vectorOfEstimatedAltitudesBelowAtmosphericInterface.push_back(
                    cartesianStateIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius_ );
    }

    // Run least squares estimation process based on selected atmosphere model
    std::vector< double > vectorOfModelSpecificParameters;
    switch ( selectedOnboardAtmosphereModel_ )
    {
    case aerodynamics::exponential_atmosphere_model:
    case aerodynamics::three_wave_atmosphere_model:
    {
        // Transform density to logarithmic space


        // Use least squares polynomial fit
//        Eigen::VectorXd estimatedAtmosphereParameters = linear_algebra::getLeastSquaresPolynomialFit( );
        break;
    }
    case aerodynamics::three_term_atmosphere_model:
    {
        break;
    }
    default:
        throw std::runtime_error( "Error in atmoshere estimation of navigation system. Atmosphere model not recognized." );
    }

    // Update atmosphere settings of onboard body map
    onboardBodyMap_.at( planetName_ )->setAtmosphereModel(
                boost::make_shared< aerodynamics::CustomConstantTemperatureAtmosphere >( selectedOnboardAtmosphereModel_,
                                                                                         215.0, 197.0, 1.3,
                                                                                         vectorOfModelSpecificParameters ) );
}

} // namespace navigation

} // namespace tudat
