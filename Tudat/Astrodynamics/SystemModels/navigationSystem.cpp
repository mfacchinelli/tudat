#include <iostream>

#include <Eigen/LU>

#include "Tudat/Astrodynamics/SystemModels/navigationSystem.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"
#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/SimulationSetup/PropagationSetup/createEnvironmentUpdater.h"

#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

namespace system_models
{

//! Function to remove errors in inertial measurement unit measurements based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromInertialMeasurementUnitMeasurement( const Eigen::Vector3d& currentInertialMeasurementUnitMeasurement,
                                                                    const Eigen::Vector3d& inertialMeasurementUnitErrors )
{
    return currentInertialMeasurementUnitMeasurement - inertialMeasurementUnitErrors;
}

//! Function to create navigation objects for onboard state estimation.
void NavigationSystem::createNavigationSystemObjects( const boost::function< Eigen::Vector3d( ) >& accelerometerMeasurementFunction )
{
    // Create function to return the current estimated rotation from the altimeter to the inertial frame
    boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAnglesCalculator =
            onboardBodyMap_.at( spacecraftName_ )->getFlightConditions( )->getAerodynamicAngleCalculator( );
    boost::function< Eigen::Quaterniond( ) > rotationFromAltimeterToInertialFrameFunction =
            boost::bind( &reference_frames::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                         aerodynamicAnglesCalculator, altimeterFrame_, reference_frames::inertial_frame );

    // Create filter object
    switch ( navigationFilterSettings_->filteringTechnique_ )
    {
    case filters::extended_kalman_filter:
    {
        navigationFilter_ = filters::createFilter< double, double >(
                    navigationFilterSettings_,
                    boost::bind( &NavigationSystem::onboardSystemModel, this, _1, _2, accelerometerMeasurementFunction ),
                    boost::bind( &NavigationSystem::onboardMeasurementModel, this, _1, _2, rotationFromAltimeterToInertialFrameFunction ),
                    boost::bind( &NavigationSystem::onboardSystemJacobian, this, _1, _2, accelerometerMeasurementFunction ),
                    boost::lambda::constant( Eigen::Matrix9d::Identity( ) ),
                    boost::bind( &NavigationSystem::onboardMeasurementJacobian, this, _1, _2, rotationFromAltimeterToInertialFrameFunction ),
                    boost::lambda::constant( Eigen::MatrixXd::Identity( 3 * altimeterPointingDirectionInAltimeterFrame_.size( ),
                                                                        3 * altimeterPointingDirectionInAltimeterFrame_.size( ) ) ) );
        break;
    }
    case filters::unscented_kalman_filter:
    {
        navigationFilter_ = filters::createFilter< double, double >(
                    navigationFilterSettings_,
                    boost::bind( &NavigationSystem::onboardSystemModel, this, _1, _2, accelerometerMeasurementFunction ),
                    boost::bind( &NavigationSystem::onboardMeasurementModel, this, _1, _2, rotationFromAltimeterToInertialFrameFunction ) );
        break;
    }
    default:
        throw std::runtime_error( "Error in setting up navigation system. The requested filtering technique is not supported." );
    }

    // Set initial time
    initialTime_ = navigationFilter_->getInitialTime( );
    currentTime_ = initialTime_;
    currentOrbitCounter_ = 0;

    // Retrieve navigation filter step size and estimated state
    navigationRefreshStepSize_ = navigationFilter_->getFilteringStepSize( );
    atmosphericNavigationRefreshStepSize_ = navigationRefreshStepSize_ / 20.0;
    Eigen::Vector9d initialEstimatedState = navigationFilter_->getCurrentStateEstimate( );

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
                2.0 * atmosphericNavigationRefreshStepSize_ / currentTime_, 25 );

    // Create object of dependent variables to save
    std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                          propagators::local_dynamic_pressure_dependent_variable, spacecraftName_, planetName_ ) );
    dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                          propagators::local_aerodynamic_heat_rate_dependent_variable, spacecraftName_, planetName_ ) );

    // Create object for propagation of spacecraft state with user-provided initial conditions
    onboardIntegratorSettings_ = boost::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
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

    // Create Jacobian interfaces
    onboardGravitationalAccelerationPartials_ = boost::make_shared< acceleration_partials::SphericalHarmonicsGravityPartial >(
                spacecraftName_, planetName_,
                boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                    onboardAccelerationModelMap_.at( spacecraftName_ ).at( planetName_ ).at( sphericalHarmonicsGravityIndex_ ) ) );
    onboardAerodynamicAccelerationPartials_ = boost::make_shared< acceleration_partials::AerodynamicAccelerationPartial >(
                boost::dynamic_pointer_cast< aerodynamics::AerodynamicAcceleration >(
                    onboardAccelerationModelMap_.at( spacecraftName_ ).at( planetName_ ).at( aerodynamicsIndex_ ) ),
                boost::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                    onboardBodyMap_.at( spacecraftName_ )->getFlightConditions( ) ),
                boost::bind( &simulation_setup::Body::getState, onboardBodyMap_.at( spacecraftName_ ) ),
                boost::bind( &simulation_setup::Body::setState, onboardBodyMap_.at( spacecraftName_ ), _1 ),
                spacecraftName_, planetName_ );
}

//! Function to post-process the accelerometer measurements.
void NavigationSystem::postProcessAccelerometerMeasurements(
        std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface )
{
    // Inform user
    std::cout << std::endl << "Removing Accelerometer Errors." << std::endl;

    // Extract current variable history
    std::vector< std::vector< double > > currentVariableHistory;
    currentVariableHistory.resize( estimatedAccelerometerErrors_.size( ) );
    std::map< double, Eigen::VectorXd > currentOrbitNavigationFilterEstimatedState = navigationFilter_->getEstimatedStateHistory( );
    for ( std::map< double, Eigen::VectorXd >::const_iterator stateHistoryIterator = currentOrbitNavigationFilterEstimatedState.begin( );
          stateHistoryIterator != currentOrbitNavigationFilterEstimatedState.end( ); stateHistoryIterator++ )
    {
        for ( unsigned int i = 0; i < currentVariableHistory.size( ); i++ )
        {
            currentVariableHistory.at( i ).push_back( stateHistoryIterator->second[ i + 6 ] );
        }
    }

    // Extract median of accelerometer errors
    for ( unsigned int i = 0; i < estimatedAccelerometerErrors_.rows( ); i++ )
    {
        estimatedAccelerometerErrors_[ i ] = statistics::computeSampleMedian( currentVariableHistory.at( i ) );
    }

    // Remove errors from accelerometer measurements and convert to inertial frame
    for ( unsigned int i = 0; i < vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
    {
        vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ) =
                removeErrorsFromInertialMeasurementUnitMeasurement(
                    vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ), estimatedAccelerometerErrors_ );
    }

    // Apply smoothing method to noisy accelerometer data
    unsigned int numberOfSamplePoints = static_cast< unsigned int >( 50.0 / atmosphericNavigationRefreshStepSize_ );
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

    input_output::writeDataMapToTextFile( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                          "kepler_" + std::to_string( currentOrbitCounter_ ) + ".dat", "/Users/Michele/Desktop/Results/"  );
    input_output::writeMatrixToFile(
                utilities::convertStlVectorToEigenVector( vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ),
                "aero_" + std::to_string( currentOrbitCounter_ ) + ".dat", 16, "/Users/Michele/Desktop/Results/" );

    // Determine actual periapse time
    double estimatedActualPeriapseTime;
    try
    {
        // Set root-finder function as the area below the acceleration curve
        estimatedActualPeriapseTime = areaBisectionRootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                        boost::bind( &areaBisectionFunction, _1, atmosphericNavigationRefreshStepSize_, timesBelowAtmosphericInterface,
                                     vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ) ) );
    }
    catch ( std::runtime_error& caughtException )
    {
        // Inform user on error
        std::cerr << "Error while computing actual periapse time. Caught this exception during root-finder "
                     "operation: " << caughtException.what( ) << std::endl
                  << "The mid-point of the time vector will be used as periapsis time." << std::endl;

        // Take the mid-point as actual periapsis time
        estimatedActualPeriapseTime =
                timesBelowAtmosphericInterface[ static_cast< unsigned int >( 0.5 * timesBelowAtmosphericInterface.rows( ) ) ];
    }

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
                mathematical_constants::PI, false, true,
                boost::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 0.001, 25 ) );

    Eigen::Vector6d estimatedCartesianStateAtNextApoapsis =
            propagateTranslationalStateWithCustomTerminationSettings( terminationSettings ).second.first.rbegin( )->second;
    Eigen::Vector6d estimatedKeplerianStateAtNextApoapsis =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                estimatedCartesianStateAtNextApoapsis, planetaryGravitationalParameter_ );
    std::cout << "Propagated: " << estimatedKeplerianStateAtNextApoapsis.transpose( ) << std::endl
              << "Initial: " << currentEstimatedKeplerianState_.transpose( ) << std::endl
              << "Previous: " << estimatedKeplerianStateAtPreviousApoapsis_.transpose( ) << std::endl;

    // Compute estimated change in velocity (i.e., Delta V) due to aerodynamic acceleration
    double estimatedChangeInVelocity = - periapseEstimatorConstant_ * numerical_quadrature::performExtendedSimpsonsQuadrature(
                atmosphericNavigationRefreshStepSize_, vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
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
    estimatedErrorInKeplerianState[ 0 ] = estimatedErrorInSemiMajorAxis;
    estimatedErrorInKeplerianState[ 1 ] = estimatedErrorInEccentricity;
    //    estimatedErrorInKeplerianState[ 5 ] = estimatedErrorInTrueAnomaly;
    historyOfEstimatedErrorsInKeplerianState_[ currentOrbitCounter_ ] = estimatedErrorInKeplerianState;

    // Compute updated estimate in Keplerian state at current time by removing the estimated error
    Eigen::Vector6d updatedCurrentKeplerianState = currentEstimatedKeplerianState_;
    updatedCurrentKeplerianState -= estimatedErrorInKeplerianState;
    // the updated state is initially defined as the one with semi-major axis and eccentricity from the time before the
    // atmospheric pass and inclination, right ascension of ascending node, argument of periapsis and true anomaly from the
    // latest estimate; then the error in semi-major axis, eccentricity and true anomaly is subtracted

    // Update navigation system state estimates
    Eigen::Matrix9d currentEstimatedCovarianceMatrix = navigationFilter_->getCurrentCovarianceEstimate( );
    currentEstimatedCovarianceMatrix.block( 0, 0, 6, 6 ).setIdentity( );
    setCurrentEstimatedKeplerianState( updatedCurrentKeplerianState );//, currentEstimatedCovarianceMatrix );
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

    input_output::writeMatrixToFile(
                utilities::convertStlVectorToEigenVector( vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ),
                "acceleration" + std::to_string( currentOrbitCounter_ ) + ".dat", 16, "/Users/Michele/Desktop/Results/" );

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

        input_output::writeMatrixToFile( estimatedAtmosphericDensitiesBelowAtmosphericInterface,
                                         "density" + std::to_string( currentOrbitCounter_ ) + ".dat", 16, "/Users/Michele/Desktop/Results/" );
        input_output::writeMatrixToFile( estimatedAltitudesBelowAtmosphericInterface,
                                         "altitude" + std::to_string( currentOrbitCounter_ ) + ".dat", 16, "/Users/Michele/Desktop/Results/" );

        // Find periapsis altitude
        double estimatedPeriapsisAltitude = estimatedAltitudesBelowAtmosphericInterface.minCoeff( );

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
                        initialParameterEstimates, estimatedAtmosphericDensitiesBelowAtmosphericInterface.array( ).log( ), 1e-6, 1e-5, 100 );

            // Add reference altitude to list of parameters and revert from logarithmic space
            modelSpecificParameters.resize( 6 );
            modelSpecificParameters[ 0 ] = estimatedPeriapsisAltitude;
            modelSpecificParameters[ 1 ] = std::exp( estimatedAtmosphereModelParameters[ 0 ] );
            modelSpecificParameters[ 2 ] = 1.0 / estimatedAtmosphereModelParameters[ 1 ];
            modelSpecificParameters.segment( 3, 3 ) = estimatedAtmosphereModelParameters.segment( 2, 3 );
            break;
        }
        }

        // Check that estimated atmospheric parameters make sense
        if ( ( modelSpecificParameters[ 1 ] < 0 ) || // density at zero altitude is negative
             ( modelSpecificParameters[ 2 ] < 0 ) ) // scale height is negative
        {
            // Inform user
            std::cerr << "Warning in atmosphere estimator. The estimated atmospheric parameters are not in line "
                         "with what is expected." << std::endl
                      << "Erroneous parameters are: " << modelSpecificParameters.transpose( ) << std::endl
                      << "Taking an average of the previous estimated values, instead." << std::endl;

            // Take average of all previous values
            modelSpecificParameters.setZero( );
            for ( std::map< unsigned int, Eigen::VectorXd >::const_iterator mapIterator = historyOfEstimatedAtmosphereParameters_.begin( );
                  mapIterator != historyOfEstimatedAtmosphereParameters_.end( ); mapIterator++ )
            {
                modelSpecificParameters += mapIterator->second;
            }
            modelSpecificParameters /= historyOfEstimatedAtmosphereParameters_.size( );
        }
        std::cout << "Atmosphere values: " << modelSpecificParameters.transpose( ) << std::endl;

        // Add values to hisotry
        historyOfEstimatedAtmosphereParameters_[ currentOrbitCounter_ ] = modelSpecificParameters;

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
            unsigned int i = 0;
            for ( std::map< unsigned int, Eigen::VectorXd >::const_reverse_iterator
                  mapIterator = historyOfEstimatedAtmosphereParameters_.rbegin( );
                  mapIterator != historyOfEstimatedAtmosphereParameters_.rend( ); mapIterator++, i++ )
            {
                if ( ( i != 0 ) && ( i < numberOfSamplesForMovingAverage ) )
                {
                    modelSpecificParameters += mapIterator->second;
                }
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

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector9d NavigationSystem::onboardSystemModel(
        const double currentTime, const Eigen::Vector9d& currentEstimatedState,
        const boost::function< Eigen::Vector3d( ) >& accelerometerMeasurementFunction )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare state derivative vector
    Eigen::Vector9d currentStateDerivative = Eigen::Vector9d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentEstimatedState.segment( 3, 3 );

    // Translational dynamics
    currentStateDerivative.segment( 3, 3 ) =
            getCurrentEstimatedGravitationalTranslationalAcceleration( currentEstimatedState.segment( 0, 6 ) ) +
            removeErrorsFromInertialMeasurementUnitMeasurement( accelerometerMeasurementFunction( ),
                                                                currentEstimatedState.segment( 6, 3 ) );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::VectorXd NavigationSystem::onboardMeasurementModel(
        const double currentTime, const Eigen::Vector9d& currentEstimatedState,
        const boost::function< Eigen::Quaterniond( ) > rotationFromAltimeterToInertialFrameFunction )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::VectorXd currentMeasurement;
    currentMeasurement.resize( 3 * altimeterPointingDirectionInAltimeterFrame_.size( ) );

    // Get current estimated state
    Eigen::Vector3d currentRadialVector = currentEstimatedState.segment( 0, 3 );
    double currentRadialDistance = currentRadialVector.norm( );

    // Pre-allocate variables
    Eigen::Quaterniond currentEstimatedRotationFromAltimeterToInertialFrame = rotationFromAltimeterToInertialFrameFunction( );
    Eigen::Vector3d altimeterPointingDirectionInInertialFrame;
    double currentDotProduct;
    double currentPseudoAltitude;

    // Loop over each altimeter
    //    std::cout << "Model altitude: " << std::endl;
    //    std::cout << currentEstimatedState.segment( 0, 6 ).transpose( ) << std::endl;
    for ( unsigned int i = 0; i < altimeterPointingDirectionInAltimeterFrame_.size( ); i++ )
    {
        // Transform pointing direction to inertial frame
        altimeterPointingDirectionInInertialFrame = currentEstimatedRotationFromAltimeterToInertialFrame *
                altimeterPointingDirectionInAltimeterFrame_.at( i );

        // Compute dot product of radial distance and pointing direction (cosine of pointing angle)
        currentDotProduct = currentRadialVector.dot( altimeterPointingDirectionInInertialFrame );

        // Determine pseudo-altitude measurement
        currentPseudoAltitude = currentDotProduct - std::sqrt(
                    planetaryRadius_ * planetaryRadius_ - currentRadialDistance * currentRadialDistance +
                    currentDotProduct * currentDotProduct );

        // Transform pseudo-altitude in vector format (in inertial frame)
        currentMeasurement.segment( 3 * i, 3 ) = currentPseudoAltitude * altimeterPointingDirectionInInertialFrame;
        //        std::cout << i << ": " << currentMeasurement.segment( 3 * i, 3 ).norm( ) << ", " <<
        //                     currentMeasurement.segment( 3 * i, 3 ).transpose( ) << std::endl
        //                  << altimeterPointingDirectionInInertialFrame.transpose( ) << std::endl;
    }

    // Return quaternion vector
    return currentMeasurement;
}

//! Function to model the onboard system Jacobian based on the simplified onboard model.
Eigen::Matrix9d NavigationSystem::onboardSystemJacobian(
        const double currentTime, const Eigen::Vector9d& currentEstimatedState,
        const boost::function< Eigen::Vector3d( ) >& accelerometerMeasurementFunction )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare Jacobian matrix and set to zero
    Eigen::Matrix9d currentSystemJacobian = Eigen::Matrix9d::Zero( );

    // Add terms due to velocity
    currentSystemJacobian( 0, 3 ) = 1.0;
    currentSystemJacobian( 1, 4 ) = 1.0;
    currentSystemJacobian( 2, 5 ) = 1.0;

    // Add terms due to spherical harmonics acceleration
    currentSystemJacobian.block( 3, 0, 3, 3 ) = onboardGravitationalAccelerationPartials_->getCurrentSphericalHarmonicsAccelerationPartial( );

    // Add terms due to aerodynamic acceleration
    currentSystemJacobian.block( 3, 0, 3, 6 ) += onboardAerodynamicAccelerationPartials_->getCurrentAerodynamicAccelerationPartial( );

    // Add terms due to accelerometer bias error
    currentSystemJacobian( 3, 6 ) = 1.0;
    currentSystemJacobian( 4, 7 ) = 1.0;
    currentSystemJacobian( 5, 8 ) = 1.0;

    // Give output
    return currentSystemJacobian;
}

//! Function to model the onboard measurements Jacobian based on the simplified onboard model.
Eigen::MatrixXd NavigationSystem::onboardMeasurementJacobian(
        const double currentTime, const Eigen::Vector9d& currentEstimatedState,
        const boost::function< Eigen::Quaterniond( ) > rotationFromAltimeterToInertialFrameFunction )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare Jacobian matrix and set to zero
    Eigen::MatrixXd currentMeasurementJacobian;
    currentMeasurementJacobian.resize( 3 * altimeterPointingDirectionInAltimeterFrame_.size( ), 12 );
    currentMeasurementJacobian.setZero( );

    // Get current estimated state
    Eigen::Vector3d currentRadialVector = currentEstimatedState.segment( 0, 3 );
    double currentRadialDistance = currentRadialVector.norm( );

    // Pre-allocate variables
    Eigen::Quaterniond currentEstimatedRotationFromAltimeterToInertialFrame = rotationFromAltimeterToInertialFrameFunction( );
    Eigen::Vector3d altimeterPointingDirectionInInertialFrame;
    double denominator;
    double currentDotProduct;

    // Loop over each altimeter
    for ( unsigned int i = 0; i < altimeterPointingDirectionInAltimeterFrame_.size( ); i++ )
    {
        // Transform pointing direction to inertial frame
        altimeterPointingDirectionInInertialFrame = currentEstimatedRotationFromAltimeterToInertialFrame *
                altimeterPointingDirectionInAltimeterFrame_.at( i );

        // Compute dot product of radial distance and pointing direction (cosine of pointing angle)
        currentDotProduct = currentRadialVector.dot( altimeterPointingDirectionInInertialFrame );
        denominator = std::sqrt( planetaryRadius_ * planetaryRadius_ -
                                 currentRadialDistance * currentRadialDistance + currentDotProduct * currentDotProduct );

        // Compute Jacobian
        double currentPseudoAltitudeDerivative;
        for ( unsigned int n = 0; n < 3; n++ )
        {
            currentPseudoAltitudeDerivative = altimeterPointingDirectionInInertialFrame[ n ] - (
                        currentDotProduct * altimeterPointingDirectionInInertialFrame[ n ] -
                        currentEstimatedState[ n ] ) / denominator;
            for ( unsigned int m = 0; m < 3; m++ )
            {
                currentMeasurementJacobian( 3 * i + m, n ) = altimeterPointingDirectionInInertialFrame[ m ] * currentPseudoAltitudeDerivative;
            }
        }
    }

    // Give output
    //    std::cout << currentMeasurementJacobian << std::endl;
    return currentMeasurementJacobian;
}

} // namespace system_models

} // namespace tudat
