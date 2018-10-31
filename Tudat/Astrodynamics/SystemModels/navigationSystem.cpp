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
    return ( currentInertialMeasurementUnitMeasurement - inertialMeasurementUnitErrors );
}

//! Function to create navigation objects for onboard state estimation.
void NavigationSystem::createNavigationSystemObjects( const unsigned int saveFrequency )
{
    // Create filter object
    switch ( navigationFilterSettings_->filteringTechnique_ )
    {
    case filters::extended_kalman_filter:
    {
        navigationFilter_ = filters::createFilter< double, double >( navigationFilterSettings_,
                                                                     boost::bind( &NavigationSystem::onboardSystemModel, this, _1, _2 ),
                                                                     boost::bind( &NavigationSystem::onboardMeasurementModel, this, _1, _2 ),
                                                                     boost::bind( &NavigationSystem::onboardSystemJacobian, this, _1, _2 ),
                                                                     boost::lambda::constant( Eigen::Matrix10d::Identity( ) ),
                                                                     boost::bind( &NavigationSystem::onboardMeasurementJacobian, this, _1, _2 ),
                                                                     boost::lambda::constant( Eigen::Matrix3d::Identity( ) ) );
        break;
    }
    case filters::unscented_kalman_filter:
    {
        navigationFilter_ = filters::createFilter< double, double >( navigationFilterSettings_,
                                                                     boost::bind( &NavigationSystem::onboardSystemModel, this, _1, _2 ),
                                                                     boost::bind( &NavigationSystem::onboardMeasurementModel, this, _1, _2 ) );
        break;
    }
    default:
        throw std::runtime_error( "Error in setting up navigation system. The requested filtering technique is not supported." );
    }

    // Set time-related parameters
    initialTime_ = navigationFilter_->getInitialTime( );
    currentTime_ = initialTime_;
    currentOrbitCounter_ = 0;
    saveFrequency_ = saveFrequency;
    saveIndex_ = 0;

    // Retrieve navigation filter step size and estimated state
    navigationRefreshStepSize_ = navigationFilter_->getFilteringStepSize( );
    Eigen::Vector10d initialNavigationFilterState = navigationFilter_->getCurrentStateEstimate( );
    
    // Set atmospheric phase step size
    double areaBisectionTimeRelativeTolerance;
    if ( navigationTesting_ )
    {
        atmosphericNavigationRefreshStepSize_ = navigationRefreshStepSize_;
        areaBisectionTimeRelativeTolerance = 2.0 * atmosphericNavigationRefreshStepSize_ / currentTime_ / 5.0;
    }
    else
    {
        atmosphericNavigationRefreshStepSize_ = navigationRefreshStepSize_ / 5.0;
        areaBisectionTimeRelativeTolerance = 2.0 * atmosphericNavigationRefreshStepSize_ / currentTime_;
    }

    // Set initial translational state
    setCurrentEstimatedCartesianState( initialNavigationFilterState.segment( cartesian_position_index, 6 ) );
    // this function also automatically stores the full state estimates at the current time

    // Store the apoapsis value of Keplerian state for the Periapse Time Estimator
    estimatedKeplerianStateAtPreviousApoapsis_ = currentEstimatedKeplerianState_;

    // Update body and acceleration maps
    updateOnboardModel( );

    // Create root-finder object for bisection of aerodynamic acceleration curve
    // The values inserted are the tolerance in independent value (i.e., about twice the difference between
    // two time steps) and the maximum number of iterations (i.e., 25 iterations)
    areaBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( areaBisectionTimeRelativeTolerance, 25 );

    // Create object of dependent variables to save
    std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                          propagators::local_dynamic_pressure_dependent_variable, spacecraftName_, planetName_ ) );
    dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                          propagators::local_aerodynamic_heat_rate_dependent_variable, spacecraftName_, planetName_ ) );

    // Create object for propagation of spacecraft state with user-provided initial conditions
    onboardIntegratorSettings_ = boost::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                0.0, 10.0, numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg56, 0.1, 100.0, 1.0e-12, 1.0e-12 );
    onboardPropagatorSettings_ = boost::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                std::vector< std::string >{ planetName_ }, onboardAccelerationModelMap_,
                std::vector< std::string >{ spacecraftName_ }, Eigen::Vector6d::Zero( ),
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

//! Function to improve the estimate of the translational state when transitioning from aided to unaided navigation.
void NavigationSystem::improveStateEstimateOnNavigationPhaseTransition( )
{
    // Inform user
    std::cout << std::endl << "Improving State Estimate for Better Phase Transition." << std::endl;

    // Declare improved state vector
    Eigen::Vector6d improvedEstimatedState;

    // Check to see if there are enough elements
    unsigned int numberOfSamplePoints = static_cast< unsigned int >( 60.0 / navigationRefreshStepSize_ );
    numberOfSamplePoints = ( ( numberOfSamplePoints % 2 ) == 0 ) ? numberOfSamplePoints + 1 : numberOfSamplePoints; // only odd values
    std::map< double, Eigen::VectorXd > historyOfEstimatedStates = navigationFilter_->getEstimatedStateHistory( );
    if ( historyOfEstimatedStates.size( ) > numberOfSamplePoints )
    {
        // Extract history of estimated states and erase extra elements
        unsigned int i = 0;
        for ( std::map< double, Eigen::VectorXd >::reverse_iterator stateHistoryIterator = historyOfEstimatedStates.rbegin( );
              stateHistoryIterator != historyOfEstimatedStates.rend( ); stateHistoryIterator++, i++ )
        {
            if ( i > numberOfSamplePoints )
            {
                historyOfEstimatedStates.erase( historyOfEstimatedStates.begin( ), stateHistoryIterator.base( ) );
                break;
            }
        }
        double timeAtBeginningOfSampling = historyOfEstimatedStates.begin( )->first;

        // Predefine variables
        double currentRelativeTime;
        Eigen::VectorXd historyOfTrueAnomalies;
        historyOfTrueAnomalies.resize( historyOfEstimatedStates.size( ) );
        Eigen::VectorXd historyOfRelativeTimes;
        historyOfRelativeTimes.resize( historyOfEstimatedStates.size( ) );

        std::vector< std::vector< double > > estimatedKeplerianStates;
        estimatedKeplerianStates.resize( 5 ); // ignore true anomaly

        // Extract information from history of states
        i = 0; // overwrite
        Eigen::Vector6d currentState;
        for ( std::map< double, Eigen::VectorXd >::const_iterator stateHistoryIterator = historyOfEstimatedStates.begin( );
              stateHistoryIterator != historyOfEstimatedStates.end( ); stateHistoryIterator++, i++ )
        {
            // Get current Keplerian state
            currentState = orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                        stateHistoryIterator->second.segment( cartesian_position_index, 6 ), planetaryGravitationalParameter_ );

            // Store true anomaly and time for computation of least squares
            historyOfTrueAnomalies[ i ] = currentState[ 5 ];
            historyOfRelativeTimes[ i ] = stateHistoryIterator->first - timeAtBeginningOfSampling;

            // Store first five elements in vector for computation of median
            for ( unsigned int j = 0; j < 5; j++ )
            {
                estimatedKeplerianStates.at( j ).push_back( currentState[ j ] );
            }
        }

        // Improve first five Keplerian elements by using their median
        for ( unsigned int i = 0; i < 5; i++ )
        {
            improvedEstimatedState[ i ] = statistics::computeSampleMedian( estimatedKeplerianStates.at( i ) );
        }

        // Use least squares to estimate the coefficient of the true anomaly approximation function
        std::vector< double > vectorOfPolynomialPowers = { 0, 1, 2 };
        Eigen::VectorXd estimatedTrueAnomalyFunctionParameters = linear_algebra::getLeastSquaresPolynomialFit(
                    historyOfRelativeTimes, historyOfTrueAnomalies, vectorOfPolynomialPowers );

        // Improve true anomaly by using the least squares estimate
        double improvedTrueAnomaly = 0.0;
        currentRelativeTime = currentTime_ - timeAtBeginningOfSampling;
        for ( unsigned int i = 0; i < vectorOfPolynomialPowers.size( ); i++ )
        {
            improvedTrueAnomaly += std::pow( currentRelativeTime, vectorOfPolynomialPowers.at( i ) ) *
                    estimatedTrueAnomalyFunctionParameters[ i ];
        }
        improvedEstimatedState[ 5 ] = currentEstimatedKeplerianState_[ 5 ];//improvedTrueAnomaly;

        // Set new value of state
        setCurrentEstimatedKeplerianState( improvedEstimatedState );
        updateOnboardModel( );
    }
    else
    {
        // Inform user
        std::cerr << "Not enought data is available. Taking latest estimate, instead." << std::endl;
    }
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
    unsigned int numberOfSamplePoints = static_cast< unsigned int >( 60.0 / atmosphericNavigationRefreshStepSize_ );
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

    // If testing, store results so that they can be compared to MATLAB
    if ( navigationTesting_ )
    {
        input_output::writeDataMapToTextFile( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface, "kepler_est.dat", "PTE&AEResults/" );
        input_output::writeMatrixToFile(
                    utilities::convertStlVectorToEigenVector( vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ),
                    "aero.dat", 16, "PTE&AEResults/" );
    }

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
    // note that this represents directly the error in estimated true anomaly, since the true anomaly at
    // periapsis is zero by definition

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

    // Compute estimated change in eccentricity due to change in velocity
    // The same assumption as for the case above holds
    double estimatedChangeInEccentricityDueToChangeInVelocity = 2.0 * std::sqrt( estimatedKeplerianStateAtPreviousApoapsis_[ 0 ] *
            ( 1.0 - std::pow( estimatedKeplerianStateAtPreviousApoapsis_[ 1 ], 2 ) ) / planetaryGravitationalParameter_ ) *
            estimatedChangeInVelocity;
    std::cout << "Estimated Change in Eccentricity: " <<
                 estimatedChangeInEccentricityDueToChangeInVelocity << std::endl;

    // Store estimated change in Keplerian elements
    Eigen::Vector6d estimatedChangeInKeplerianState = Eigen::Vector6d::Zero( );
    estimatedChangeInKeplerianState[ 0 ] = estimatedChangeInSemiMajorAxisDueToChangeInVelocity;
    estimatedChangeInKeplerianState[ 1 ] = estimatedChangeInEccentricityDueToChangeInVelocity;
    estimatedChangeInKeplerianState[ 5 ] = estimatedErrorInTrueAnomaly;
    historyOfEstimatedChangesInKeplerianState_[ currentOrbitCounter_ ] = estimatedChangeInKeplerianState;

    // Compute updated estimate in Keplerian state at current time by removing the estimated change in elements
    Eigen::Vector6d updatedCurrentKeplerianState;
    updatedCurrentKeplerianState.segment( 0, 2 ) = estimatedKeplerianStateAtPreviousApoapsis_.segment( 0, 2 );
    updatedCurrentKeplerianState.segment( 3, 4 ) = currentEstimatedKeplerianState_.segment( 3, 4 );
    updatedCurrentKeplerianState += estimatedChangeInKeplerianState;
    // the updated state is initially defined as the one with semi-major axis and eccentricity from the apoapsis before the
    // atmospheric pass and inclination, right ascension of ascending node, argument of periapsis and true anomaly from the
    // latest estimate; then the change in semi-major axis, eccentricity and true anomaly is added

    // Update navigation system state estimates
    setCurrentEstimatedKeplerianState( updatedCurrentKeplerianState );

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

        // If testing, store results so that they can be compared to MATLAB
        if ( navigationTesting_ )
        {
            input_output::writeMatrixToFile(
                        utilities::convertStlVectorToEigenVector( vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ),
                        "acceleration.dat", 16, "PTE&AEResults/" );
            input_output::writeMatrixToFile( estimatedAtmosphericDensitiesBelowAtmosphericInterface, "density.dat", 16, "PTE&AEResults/" );
            input_output::writeMatrixToFile( estimatedAltitudesBelowAtmosphericInterface, "altitude.dat", 16, "PTE&AEResults/" );
        }

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

        // Add values to hisotry
        historyOfEstimatedAtmosphereParameters_[ currentOrbitCounter_ ] = modelSpecificParameters;

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
    else
    {
        // Inform user
        std::cerr << "Warning in atmosphere estimator. No altitudes recorded below the atmospheric interface. "
                  << "Reduced atmospheric interface: " << ( reducedAtmosphericInterfaceRadius_ - planetaryRadius_ ) / 1.0e3 << " km. "
                  << "Lowest recorded altitude: " << utilities::convertStlVectorToEigenVector(
                         vectorOfEstimatedAltitudesBelowAtmosphericInterface ).minCoeff( ) / 1.0e3 << " km." << std::endl;
    }
}

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector10d NavigationSystem::onboardSystemModel(
        const double currentTime, const Eigen::Vector10d& currentNavigationFilterState )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare state derivative vector
    Eigen::Vector10d currentStateDerivative = Eigen::Vector10d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentNavigationFilterState.segment( cartesian_velocity_index, 3 );

    // Translational dynamics
    currentStateDerivative.segment( 3, 3 ) = getCurrentEstimatedTranslationalAcceleration(
                currentNavigationFilterState.segment( cartesian_position_index, 6 ), currentNavigationFilterState[ drag_coefficient_index ] );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector3d NavigationSystem::onboardMeasurementModel(
        const double currentTime, const Eigen::Vector10d& currentNavigationFilterState )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::Vector3d currentMeasurementVector;

    // Add terms due to accelerometer bias error
    currentMeasurementVector = currentNavigationFilterState.segment( 6, 3 );

    // Add terms due to aerodynamic acceleration
    currentMeasurementVector += getCurrentEstimatedNonGravitationalTranslationalAcceleration(
                currentNavigationFilterState.segment( cartesian_position_index, 6 ), currentNavigationFilterState[ drag_coefficient_index ] );

    // Give output
    return currentMeasurementVector;
}

//! Function to model the onboard system Jacobian based on the simplified onboard model.
Eigen::Matrix10d NavigationSystem::onboardSystemJacobian(
        const double currentTime, const Eigen::Vector10d& currentNavigationFilterState )
{
    TUDAT_UNUSED_PARAMETER( currentTime );
    updateOnboardModel( currentNavigationFilterState.segment( cartesian_position_index, 6 ), currentNavigationFilterState[ drag_coefficient_index ] );

    // Declare Jacobian matrix and set to zero
    Eigen::Matrix10d currentSystemJacobian = Eigen::Matrix10d::Zero( );

    // Add terms due to velocity
    currentSystemJacobian( 0, 3 ) = 1.0;
    currentSystemJacobian( 1, 4 ) = 1.0;
    currentSystemJacobian( 2, 5 ) = 1.0;

    // Add terms due to spherical harmonics acceleration
    currentSystemJacobian.block( 3, 0, 3, 3 ) +=
            onboardGravitationalAccelerationPartials_->getCurrentSphericalHarmonicsAccelerationPartial( );

    // Add terms due to aerodynamic acceleration
    currentSystemJacobian.block( 3, 0, 3, 6 ) += onboardAerodynamicAccelerationPartials_->getCurrentAerodynamicAccelerationPartial( );
    currentSystemJacobian.block( 3, 9, 3, 1 ) +=
            onboardAerodynamicAccelerationPartials_->getCurrentAerodynamicAccelerationPartialWrtDragCoefficient( );

    // Give output
    return currentSystemJacobian;
}

//! Function to model the onboard measurements Jacobian based on the simplified onboard model.
Eigen::Matrix< double, 3, 10 > NavigationSystem::onboardMeasurementJacobian(
        const double currentTime, const Eigen::Vector10d& currentNavigationFilterState )
{
    TUDAT_UNUSED_PARAMETER( currentTime );
    updateOnboardModel( currentNavigationFilterState.segment( cartesian_position_index, 6 ), currentNavigationFilterState[ drag_coefficient_index ] );

    // Declare Jacobian matrix and set to zero
    Eigen::Matrix< double, 3, 10 > currentMeasurementJacobian = Eigen::Matrix< double, 3, 10 >::Zero( );

    // Add terms due to aerodynamic acceleration
    currentMeasurementJacobian.block( 0, 0, 3, 6 ) += onboardAerodynamicAccelerationPartials_->getCurrentAerodynamicAccelerationPartial( );
    currentMeasurementJacobian.block( 0, 9, 3, 1 ) +=
            onboardAerodynamicAccelerationPartials_->getCurrentAerodynamicAccelerationPartialWrtDragCoefficient( );

    // Add terms due to accelerometer bias error
    currentMeasurementJacobian( 0, 6 ) = 1.0;
    currentMeasurementJacobian( 1, 7 ) = 1.0;
    currentMeasurementJacobian( 2, 8 ) = 1.0;

    // Give output
    return currentMeasurementJacobian;
}

} // namespace system_models

} // namespace tudat
