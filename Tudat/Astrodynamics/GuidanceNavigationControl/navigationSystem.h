/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NAVIGATION_SYSTEM_H
#define TUDAT_NAVIGATION_SYSTEM_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/customConstantTemperatureAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/extraFunctions.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"
#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"

#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Mathematics/Filters/createFilter.h"

//! Typedefs and using statements to simplify code.
namespace Eigen { typedef Eigen::Matrix< double, 12, 1 > Vector12d; typedef Eigen::Matrix< double, 5, 5 > Matrix5d;
                  typedef Eigen::Matrix< double, 12, 12 > Matrix12d; }

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to remove errors in inertial measurement unit measurements based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromInertialMeasurementUnitMeasurement( const Eigen::Vector3d& currentInertialMeasurementUnitMeasurement,
                                                                    const Eigen::Vector6d& inertialMeasurementUnitErrors );

//! Function to remove errors in the altimeter measurements based on the current orientation.
void removeErrorsFromAltimeterMeasurement( double& currentMeasuredAltitude, const double planetaryRadius );

//! Class for the navigation system of an aerobraking maneuver.
class NavigationSystem
{
public:

    //! Enumeration for indices of navigation state vector.
    enum NavigationStateIndices
    {
        cartesian_position_index = 0,
        cartesian_velocity_index = 3,
        quaternion_real_index = 6,
        quaternion_imaginary_index = 7,
        accelerometer_bias_index = 10,
        accelerometer_scale_factor_index = 13,
        gyroscope_bias_index = 16,
        gyroscope_scale_factor_index = 19
    };

    //! Enumeration for navigation phases.
    enum NavigationPhaseIndicator
    {
        unaided_navigation_phase = 0,
        optical_navigation_phase = 1,
        altimeter_navigation_phase = 2
    };

    //! Constructor.
    /*!
     *  Constructor for the navigation system of an aerobraking maneuver.
     *  \param navigationFilter
     *  \param stateTransitionMatrixFunction Function to compute the state transition matrix at the current time and state.
     */
    NavigationSystem( const simulation_setup::NamedBodyMap& onboardBodyMap,
                      const basic_astrodynamics::AccelerationMap& onboardAccelerationModelMap,
                      const std::string& spacecraftName, const std::string& planetName,
                      const boost::shared_ptr< filters::FilterSettings< > > navigationFilterSettings,
                      const aerodynamics::AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel,
                      const double atmosphericInterfaceAltitude, const double reducedAtmosphericInterfaceAltitude,
                      const unsigned int numberOfRequiredAtmosphereSamplesForInitiation,
                      const Eigen::Vector3d& altimeterBodyFixedPointingDirection,
                      const std::pair< double, double > altimeterAltitudeRange ) :
        onboardBodyMap_( onboardBodyMap ), onboardAccelerationModelMap_( onboardAccelerationModelMap ),
        spacecraftName_( spacecraftName ), planetName_( planetName ), navigationFilterSettings_( navigationFilterSettings ),
        selectedOnboardAtmosphereModel_( selectedOnboardAtmosphereModel ),
        planetaryGravitationalParameter_( onboardBodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( ) ),
        planetaryRadius_( onboardBodyMap_.at( planetName_ )->getShapeModel( )->getAverageRadius( ) ),
        atmosphericInterfaceRadius_( planetaryRadius_ + atmosphericInterfaceAltitude ),
        reducedAtmosphericInterfaceRadius_( planetaryRadius_ + reducedAtmosphericInterfaceAltitude ),
        numberOfRequiredAtmosphereSamplesForInitiation_( numberOfRequiredAtmosphereSamplesForInitiation ),
        altimeterBodyFixedPointingDirection_( altimeterBodyFixedPointingDirection ), altimeterAltitudeRange_( altimeterAltitudeRange )
    {
        // Create environment updater
        createOnboardEnvironmentUpdater( );

        // Atmosphere estimator initialization
        atmosphereEstimatorInitialized_ = false;

        // Set initial estimated accelerometer errors to zero
        estimatedAccelerometerErrors_.setZero( );

        // Get index of central body acceleration (which is not measured by the IMUs)
        for ( accelerationMapIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                if ( ( basic_astrodynamics::getAccelerationModelType( accelerationMapIterator_->second[ i ] ) ==
                       basic_astrodynamics::spherical_harmonic_gravity ) &&
                     ( accelerationMapIterator_->first == planetName_ ) )
                {
                    sphericalHarmonicsGravityIndex_ = i;
                    break;
                }
            }
        }
    }

    //! Destructor.
    ~NavigationSystem( ) { }

    //! Function to create navigation objects for onboard state estimation.
    /*!
     *  Function to create navigation objects for onboard state estimation. This function should be called before any feature of
     *  the navigation system is used, as it creates most of the objects that are needed for state estimation (i.e., navigation
     *  filter, root-finder, onboard integrator and propagator settings).
     *  \param onboardSystemModel Function that defines the onboard system model (i.e., Cowell equations of motion).
     *  \param onboardMeasurementModel Function that defines the onboard measurement model (i.e., simulated IMU and star tracker
     *      measurements).
     */
    void createNavigationSystemObjects(
            const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardSystemModel,
            const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardMeasurementModel );

    //! Function to determine the navigation phase.
    void determineNavigationPhase( )
    {
        // Declare navigation phase indicator and set value to unaided phase
        NavigationPhaseIndicator detectedNavigationPhase = unaided_navigation_phase;

        // Store previous navigation phase
        previousNavigationPhase_ = currentNavigationPhase_;

        // Set current navigation phase
        currentNavigationPhase_ = detectedNavigationPhase;
    }

    //! Function to run the State Estimator (SE).
    void runStateEstimator( const double currentTime, Eigen::Vector3d& currentExternalMeasurementVector )
    {
        // Set time
        currentTime_ = currentTime;

        // Update filter
        navigationFilter_->updateFilter( currentTime_, currentExternalMeasurementVector );

        // Extract estimated state and update navigation estimates
        Eigen::Vector12d updatedEstimatedState = navigationFilter_->getCurrentStateEstimate( );
        setCurrentEstimatedCartesianState( updatedEstimatedState.segment( 0, 6 ) );
        std::cout << "ONB Pos: " << currentEstimatedCartesianState_.transpose( ) << std::endl;

        // Update body and acceleration maps
        updateOnboardModel( );
    }

    //! Function to run post-atmosphere processes.
    /*!
     *  Function to run post-atmosphere processes. This function takes the map of full accelerations measured by the IMU and
     *  after processing them, it feeds them to the functions that run the post-atosphere processes. The processing consists of
     *  reducing the accelerations to only the values that are measured below the atmospheric interface altitude. The
     *  post-atmosphere processes are carried out in the following order:
     *      - post-process accelerations: calibrate (if atmosphere estimator is not initialized) and remove accelerometer errors
     *      - run PTE only if atmosphere estimator (AE) is initialized
     *      - run AE to estimate atmospheric properties
     *  \param mapOfMeasuredAerodynamicAcceleration Map of time and measured aerodynamic acceleration as measured by the IMU, i.e.,
     *      before any errors have been removed.
     */
    void runPostAtmosphereProcesses( const std::map< double, Eigen::Vector3d >& mapOfMeasuredAerodynamicAcceleration )
    {
        using mathematical_constants::PI;

        // Extract aerodynamic accelerations of when the spacecraft is below the atmospheric interface altitude
        double currentIterationTime;
        std::map< double, Eigen::Vector6d > mapOfEstimatedCartesianStatesBelowAtmosphericInterface;
        std::map< double, Eigen::Vector6d > mapOfEstimatedKeplerianStatesBelowAtmosphericInterface;
        std::map< double, Eigen::Vector3d > mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface;
        std::vector< Eigen::Vector3d > vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface;
        for ( translationalStateIterator_ = currentOrbitHistoryOfEstimatedTranslationalStates_.begin( );
              translationalStateIterator_ != currentOrbitHistoryOfEstimatedTranslationalStates_.end( ); translationalStateIterator_++ )
        {
            if ( translationalStateIterator_->second.first.segment( 0, 3 ).norm( ) <= atmosphericInterfaceRadius_ )
            {
                // Retireve time, state and acceleration of where the altitude is below the atmospheric interface
                currentIterationTime = translationalStateIterator_->first;
                mapOfEstimatedCartesianStatesBelowAtmosphericInterface[ currentIterationTime ] =
                        translationalStateIterator_->second.first;
                mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ] =
                        translationalStateIterator_->second.second;
                mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface[ currentIterationTime ] =
                        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_[ currentIterationTime ];
                vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.push_back(
                            mapOfMeasuredAerodynamicAcceleration.at( currentIterationTime ) );

                // Modify the true anomaly such that it is negative where it is above PI radians (before estimated periapsis)
                if ( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ][ 5 ] >= PI )
                {
                    mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ][ 5 ] -= 2.0 * PI;
                }
            }
        }

        // Only proceed if satellite flew below atmospheric interface altitude
        if ( vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ) > 0 )
        {
            // Remove errors from measured accelerations
            postProcessAccelerometerMeasurements( vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface,
                                                  mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface );
            // this function also calibrates the accelerometers if it is the first orbit

            // Compute magnitude of aerodynamic acceleration
            std::vector< double > vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface;
            for ( unsigned int i = 0; i < vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
            {
                vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back(
                            vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ).norm( ) );
            }

            // Run periapse time estimator if not the first orbit
//            if ( atmosphereEstimatorInitialized_ ) // historyOfEstimatedAtmosphereParameters_.size( ) > 0 ) //
            {
                runPeriapseTimeEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                          vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
            }

            // Run atmosphere estimator with processed results
            runAtmosphereEstimator( mapOfEstimatedCartesianStatesBelowAtmosphericInterface,
                                    vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
        }
    }

    //! Function to update current translational Cartesian state with Deep Space Network measurement.
    /*!
     *  Function to update current translational Cartesian state with Deep Space Network measurement. Due to light-time delays,
     *  the Deep Space Network (DSN) measurement is taken at an earlier time. Thus, this function propagates the received state
     *  information to the current time.
     *  \param currentDeepSpaceNetworkTrackingData Pair of double and vector, denoting the light-time delay at the time of measurement
     *      (2 times the light-time distance to the spacecraft), and the Cartesian state computed by merging the tracking result with
     *      a post-processing software for orbit determination (done on the ground).
     */
    void processDeepSpaceNetworkTracking( const std::pair< double, Eigen::Vector6d >& currentDeepSpaceNetworkTrackingData )
    {
        // Split Deep Space Network tracking data in light-time delay and Cartesian state
        double currentLightTimeDelay = currentDeepSpaceNetworkTrackingData.first;
        Eigen::Vector6d currentTrackedState = currentDeepSpaceNetworkTrackingData.second;

        // Propagate tracked state to current time
        boost::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings =
                boost::make_shared< propagators::PropagationTimeTerminationSettings >( currentTime_ );
        Eigen::Vector6d propagatedStateBasedOnTracking = propagateStateWithCustomTerminationSettings(
                    terminationSettings, currentTrackedState, currentTime_ - currentLightTimeDelay ).first.rbegin( )->second;

        // Reset navigation filter (including covariance)
        setCurrentEstimatedCartesianState( propagatedStateBasedOnTracking, Eigen::Matrix12d::Identity( ) );
    }

    //! Function to propagate translational Cartesian state to specified termination settings.
    /*!
     *  Function to propagate translational Cartesian state to specified termination settings. Note that this function does not
     *  continuously estimate the spacecraft state, but instead it is called by the navigation system or any other onboard system
     *  to propagate the current estimated state to some termination condition. The initial time and state are automatically set to
     *  the current time and current estimated translational Cartesian state, respectively, but they can be overwritten by setting
     *  them from the input.
     *  \param propagationTerminationSettings Termination settings to be used to stop the propagation.
     *  \return Pair of dynamics simulator results. The first element is the map of time and propagated state, whereas the second
     *      one is the map of time and dependent variable values.
     */
    std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > >
    propagateStateWithCustomTerminationSettings(
            const boost::shared_ptr< propagators::PropagationTerminationSettings > propagationTerminationSettings,
            const Eigen::Vector6d& initialTranslationalCartesianState = Eigen::Vector6d::Zero( ), const double initialTime = -1.0 )
    {
        // Set initial time
        if ( static_cast< int >( initialTime ) == -1 )
        {
            onboardIntegratorSettings_->initialTime_ = currentTime_;
        }
        else
        {
            onboardIntegratorSettings_->initialTime_ = initialTime;
        }

        // Set initial state
        if ( initialTranslationalCartesianState.norm( ) == 0.0 )
        {
            onboardPropagatorSettings_->resetInitialStates( currentEstimatedCartesianState_ );
        }
        else
        {
            onboardPropagatorSettings_->resetInitialStates( initialTranslationalCartesianState );
        }

        // Set propagation settings
        onboardPropagatorSettings_->resetTerminationSettings( propagationTerminationSettings );

        // Create dynamics simulator and propagate
        propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                    onboardBodyMap_, onboardIntegratorSettings_, onboardPropagatorSettings_ );

        // Retrieve results from onboard computer and systems
        std::map< double, Eigen::VectorXd > translationalStateResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariablesResult = dynamicsSimulator.getDependentVariableHistory( );

        // Give output
        return std::make_pair( translationalStateResult, dependentVariablesResult );
    }

    //! Function to retireve current time.
    double getCurrentTime( ) { return currentTime_; }

    //! Function to retireve current state.
    Eigen::Vector12d getCurrentEstimatedState( ) { return navigationFilter_->getCurrentStateEstimate( ); }

    //! Function to retrieve the history of estimated states directly from the navigation filter.
    std::map< double, Eigen::VectorXd > getHistoryOfEstimatedStatesFromNavigationFilter( )
    {
        return navigationFilter_->getEstimatedStateHistory( );
    }

    //! Function to retrieve the history of estimated covariance directly from the navigation filter.
    std::map< double, Eigen::MatrixXd > getHistoryOfEstimatedCovarianceFromNavigationFilter( )
    {
        return navigationFilter_->getEstimatedCovarianceHistory( );
    }

    //! Function to retireve current state.
    std::pair< Eigen::Vector6d, Eigen::Vector6d > getCurrentEstimatedTranslationalState( )
    {
        return std::make_pair( currentEstimatedCartesianState_, currentEstimatedKeplerianState_ );
    }

    //! Function to compute the current estimated mean motion.
    double getCurrentEstimatedMeanMotion( )
    {
        return basic_astrodynamics::computeKeplerMeanMotion(
                    currentEstimatedKeplerianState_[ 0 ], planetaryGravitationalParameter_ );
    }

    //! Function to retrieve the density at the input conditions according to the onboard model.
    double getDensityAtSpecifiedConditions( double altitude, double longitude = 0.0 )
    {
        return onboardBodyMap_.at( planetName_ )->getAtmosphereModel( )->getDensity( altitude, longitude, 0.0, 0.0 );
    }

    //! Function to retrieve current estimated translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated translational accelerations exerted on the spacecraft. The acceleration
     *  is computed by using the onboard accelerations model map.
     *  \return Vector denoting the full acceleration vector experienced by the spacecraft.
     */
    Eigen::Vector3d getCurrentEstimatedTranslationalAcceleration( )
    {
        return currentEstimatedTranslationalAcceleration_;
    }

    //! Function to retrieve current estimated gravitational translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated gravitational translational accelerations exerted on the spacecraft. The
     *  acceleration is computed by using the onboard accelerations model map. Note that this vector represents only the
     *  gravitational accelerations.
     *  \return Vector denoting only the gravitational accelerations experienced by the spacecraft.
     */
    Eigen::Vector3d getCurrentEstimatedGravitationalTranslationalAcceleration( )
    {
        return currentEstimatedGravitationalTranslationalAcceleration_;
    }

    //! Function to retrieve current estimated non-gravitational translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated non-gravitational translational accelerations exerted on the spacecraft. The
     *  acceleration is computed by using the onboard accelerations model map. Note that this vector represents only the
     *  non-gravitational accelerations (which are supposed to emulate the accelerations measured by the IMU).
     *  \return Vector denoting only the non-gravitational (i.e., aerodynamic) accelerations experienced by the spacecraft.
     */
    Eigen::Vector3d getCurrentEstimatedNonGravitationalTranslationalAcceleration( )
    {
        return currentEstimatedNonGravitationalTranslationalAcceleration_;
    }

    //! Function to retrieve the estimated accelerometer errors.
    Eigen::Vector6d getEstimatedAccelerometerErrors( )
    {
        return estimatedAccelerometerErrors_;
    }

    //! Function to retrieve information on the initialization of the atmosphere estimator.
    /*!
     *  Function to retrieve information on the initialization of the atmosphere estimator.
     *  \return Pair of integers, where the first element denotes the number of atmosphere samples that have been taken so far,
     *      and the second element indicates the number of atmosphere samples required for the atmosphere estimator to be
     *      considered initialized.
     */
    std::pair< unsigned int, unsigned int > getAtmosphereInitiationIndicators( )
    {
        return std::make_pair( historyOfEstimatedAtmosphereParameters_.size( ), numberOfRequiredAtmosphereSamplesForInitiation_ );
    }

    //! Function to retrieve history of estimated translational and rotational states over time.
    /*!
     *  Function to retrieve history of estimated translational and rotational states over time.
     *  \return Map where the keys are times and the mapped value is a pair, whose first element is the pair of
     *      Cartesian and Keplerian elements, whereas the second element is the rotational state.
     */
    std::map< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getHistoryOfEstimatedStates( )
    {
        return historyOfEstimatedStates_;
    }

    //! Function to retrieve history of estimated translational accelerations for the current orbit.
    std::map< double, Eigen::Vector3d > getCurrentOrbitHistoryOfEstimatedTranslationalAccelerations( )
    {
        return currentOrbitHistoryOfEstimatedTranslationalAccelerations_;
    }

    //! Function to retrieve history of estimated non-gravitational translational accelerations for the current orbit.
    std::map< double, Eigen::Vector3d > getCurrentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations( )
    {
        return currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_;
    }

    //! Function to retireve refresh step size of navigation system.
    double getNavigationRefreshStepSize( ) { return navigationRefreshStepSize_; }

    //! Function to retireve the standard gravitational parameter of the body being orbited.
    double getStandardGravitationalParameter( ) { return planetaryGravitationalParameter_; }

    //! Function to retireve the radius of the body being orbited.
    double getRadius( ) { return planetaryRadius_; }

    //! Function to retireve the atmospheric interface radius of the body being orbited.
    double getAtmosphericInterfaceRadius( ) { return atmosphericInterfaceRadius_; }

    //! Function to set current Cartesian state to new value.
    /*!
     *  Function to set current Cartesian state to new value. The value of the Keplerian state is set
     *  automatically by converting the Cartesian state to Keplerian elements. Also, the function updates the
     *  history of translational to the new values.
     *  \param newCurrentKeplerianState New Cartesian state at current time.
     *  \param newCurrentCovarianceMatrix New covariance matrix at current time.
     */
    void setCurrentEstimatedCartesianState( const Eigen::Vector6d& newCurrentCartesianState,
                                            const Eigen::Matrix12d& newCurrentCovarianceMatrix = Eigen::Matrix12d::Zero( ) )
    {
        // Update navigation system current value
        currentEstimatedCartesianState_ = newCurrentCartesianState;
        currentEstimatedKeplerianState_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                    newCurrentCartesianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values

        // Update navigation filter current value
        Eigen::Vector12d updatedCurrentEstimatedState = navigationFilter_->getCurrentStateEstimate( );
        updatedCurrentEstimatedState.segment( 0, 6 ) = currentEstimatedCartesianState_;
        navigationFilter_->modifyCurrentStateAndCovarianceEstimates( updatedCurrentEstimatedState,
                                                                     newCurrentCovarianceMatrix );
    }

    //! Function to set current Keplerian state to new value.
    /*!
     *  Function to set current Keplerian state to new value. The value of the Cartesian state is set
     *  automatically by converting the Keplerian state to Cartesian elements. Also, the function updates the
     *  history of translational to the new values.
     *  \param newCurrentKeplerianState New Keplerian state at current time.
     *  \param newCurrentCovarianceMatrix New covariance matrix at current time.
     */
    void setCurrentEstimatedKeplerianState( const Eigen::Vector6d& newCurrentKeplerianState,
                                            const Eigen::Matrix12d& newCurrentCovarianceMatrix = Eigen::Matrix12d::Zero( ) )
    {
        // Update navigation system current value
        currentEstimatedKeplerianState_ = newCurrentKeplerianState;
        currentEstimatedCartesianState_ = orbital_element_conversions::convertKeplerianToCartesianElements(
                    newCurrentKeplerianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values

        // Update navigation filter current value
        Eigen::Vector12d updatedCurrentEstimatedState = navigationFilter_->getCurrentStateEstimate( );
        updatedCurrentEstimatedState.segment( 0, 6 ) = currentEstimatedCartesianState_;
        navigationFilter_->modifyCurrentStateAndCovarianceEstimates( updatedCurrentEstimatedState,
                                                                     newCurrentCovarianceMatrix );
    }

    //! Clear history of estimated states and accelerations for the current orbit.
    void clearCurrentOrbitEstimationHistory( )
    {
        currentOrbitHistoryOfEstimatedTranslationalStates_.clear( );
        currentOrbitHistoryOfEstimatedTranslationalAccelerations_.clear( );
        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_.clear( );
        navigationFilter_->clearFilterHistory( );
    }

    //! Integer denoting the current orbit counter.
    unsigned int currentOrbitCounter_;

private:

    Eigen::Vector6d estimatedApoapsisKeplerianState_;

    //! Function to create the onboard environment updater.
    void createOnboardEnvironmentUpdater( );

    //! Function to update the body and acceleration map with the current time and state information.
    /*!
     *  Function to update the body and acceleration map with the current time and state information.
     */
    void updateOnboardModel( )
    {
        // Update environment
        mapOfStatesToUpdate_[ propagators::translational_state ] = currentEstimatedCartesianState_;
        onboardEnvironmentUpdater_->updateEnvironment( currentTime_, mapOfStatesToUpdate_ );

        // Loop over bodies exerting accelerations on spacecraft
        for ( accelerationMapIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                // Update acceleration model
                accelerationMapIterator_->second[ i ]->resetTime( TUDAT_NAN ); // force update by resetting time
                accelerationMapIterator_->second[ i ]->updateMembers( currentTime_ );
            }
        }

        // Update accelerations based on new translational state estimate
        updateCurrentEstimatedAccelerations( );
    }

    //! Function to update current estimated accelerations exerted on the spacecraft.
    /*!
     *  Function to update current estimated accelerations exerted on the spacecraft. The acceleration is computed by using
     *  the onboard accelerations model map, and is stored in three different vectors: one with the full acceleration, one with
     *  only gravitational accelerations, and the last one with only non-gravitational accelerations. These can be retrieved with
     *  their respective getCurrentEstimatedTranslationalXXXAcceleration function.
     */
    void updateCurrentEstimatedAccelerations( )
    {
        // Define output and set accelerations to zero
        currentEstimatedTranslationalAcceleration_.setZero( );
        currentEstimatedGravitationalTranslationalAcceleration_.setZero( );
        currentEstimatedNonGravitationalTranslationalAcceleration_.setZero( );

        // Iterate over all accelerations acting on body
        Eigen::Vector3d currentAcceleration;
        for ( accelerationMapIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each accelerations
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                // Calculate acceleration and add to state derivative
                currentAcceleration = accelerationMapIterator_->second[ i ]->getAcceleration( );
                currentEstimatedTranslationalAcceleration_ += currentAcceleration;

                // Only add the gravitational accelerations
                if ( ( i == sphericalHarmonicsGravityIndex_ ) && ( accelerationMapIterator_->first == planetName_ ) )
                {
                    // Calculate acceleration and add to state derivative
                    currentEstimatedGravitationalTranslationalAcceleration_ += currentAcceleration;
                }

                // Disregard the central gravitational accelerations, since IMUs do not measure them
                if ( !( ( i == sphericalHarmonicsGravityIndex_ ) && ( accelerationMapIterator_->first == planetName_ ) ) )
                {
                    // Calculate acceleration and add to state derivative
                    currentEstimatedNonGravitationalTranslationalAcceleration_ += currentAcceleration;
                }
            }
        }
        std::cout << std::setprecision( 20 )
                  << "ONB Acc: " << currentEstimatedTranslationalAcceleration_.transpose( ) << std::endl << std::endl;

        // Store acceleration value
        currentOrbitHistoryOfEstimatedTranslationalAccelerations_[ currentTime_ ] =
                currentEstimatedTranslationalAcceleration_;
        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_[ currentTime_ ] =
                currentEstimatedNonGravitationalTranslationalAcceleration_;
    }

    //! Function to post-process the accelerometer measurements.
    /*!
     *  Function to post-process the accelerometer measurements. This function first applies a simple smoothing method
     *  (moving average) to the filter, then estimates the errors in the accelerometer (bias and scale factor) and then
     *  transforms the results from body-fixed frame (the frame of the accelerometer) to inertial (the frame in which
     *  accelerations are needed for propagation).
     *  \param vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface Vector of estimated aerodynamic acceleration
     *      below the atmospheric interface altitude and corrupted by accelerometer errors. The acceleration is thus taken
     *      directly from the IMU.
     *  \param mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface Map of time and expected aerodynamic acceleration
     *      below the atmospheric interface altitude.
     */
    void postProcessAccelerometerMeasurements(
            std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface,
            const std::map< double, Eigen::Vector3d >& mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface );

    //! Function to run the Periapse Time Estimator (PTE).
    /*!
     *  Function to estimate the time of periapsis, by finding the centroid of the aerodynamic acceleration curve,
     *  via the bisection root-finder. Then it computes the change in velocity \f$ \Delta V \f$ by integrating the aerodynamic
     *  acceleration, and the change in semi-major axis, by assuming that the change in velocity occurs as an impulsive shot at
     *  pericenter. Then the values of \f$ \Delta \vartheta \f$ (which is derived from the change in time of periapsis crossing)
     *  and \f$ \Delta a \f$ are used to correct the current state estimate, by using the state transition matrix of the
     *  system.
     *  \param mapOfEstimatedKeplerianStatesBelowAtmosphericInterface Map of time and estimated Keplerian elements below
     *      the atmospheric interface altitude.
     *  \param mapOfMeasuredAerodynamicAcceleration Map of time and estimated aerodynamic acceleration. The acceleration is
     *      the one computed from the IMU measurements and is stored as a three-dimensional vector.
     */
    void runPeriapseTimeEstimator(
            const std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
            const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

    //! Function to run the Atmosphere Estimator (AE).
    /*!
     *  Function to estimate the atmospheric properties based on the selected onboard atmosphere model.
     *  \param mapOfEstimatedCartesianStatesBelowAtmosphericInterface Map of time and estimated Cartesian elements below
     *      the atmospheric interface altitude.
     *  \param mapOfMeasuredAerodynamicAcceleration Map of time and estimated aerodynamic acceleration. The acceleration is
     *      the one computed from the IMU measurements and is stored as a three-dimensional vector.
     */
    void runAtmosphereEstimator(
            const std::map< double, Eigen::Vector6d >& mapOfEstimatedCartesianStatesBelowAtmosphericInterface,
            const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

    //! Function to store current time and current state estimates.
    void storeCurrentTimeAndStateEstimates( )
    {
        // Store translational values for current orbit
        currentOrbitHistoryOfEstimatedTranslationalStates_[ currentTime_ ] = std::make_pair( currentEstimatedCartesianState_,
                                                                                             currentEstimatedKeplerianState_ );
        historyOfEstimatedStates_[ currentTime_ ] = currentOrbitHistoryOfEstimatedTranslationalStates_[ currentTime_ ];
    }

    //! Body map of the onboard simulation.
    const simulation_setup::NamedBodyMap onboardBodyMap_;

    //! Pointer to accelerations exerted on the spacecraft according to the onboard model.
    const basic_astrodynamics::AccelerationMap onboardAccelerationModelMap_;

    //! String denoting the name of the spacecraft body.
    const std::string spacecraftName_;

    //! String denoting the name of the planet being orbited body.
    const std::string planetName_;

    //! Pointer to the filter settings to be used to create the navigation filter for state estimation.
    const boost::shared_ptr< filters::FilterSettings< > > navigationFilterSettings_;

    //! Enumeration denoting atmosphere model being used for onboard applications.
    const aerodynamics::AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel_;

    //! Standard gravitational parameter of body being orbited.
    const double planetaryGravitationalParameter_;

    //! Radius of body being orbited.
    const double planetaryRadius_;

    //! Double denoting the atmospheric interface altitude.
    /*!
     *  Double denoting the atmospheric interface altitude. This altitude corresponds to the altitude where the aerodynamic
     *  accelerations have the same magnitude as the solar radiation pressure accelerations. At these conditions, the bias of
     *  the accelerometer is easier to quantify.
     */
    const double atmosphericInterfaceRadius_;

    //! Double denoting the reduced atmospheric interface altitude.
    /*!
     *  Double denoting the reduced atmospheric interface altitude. This altitude is used as atltitude threshold for the
     *  atmosphere estimation. The reason why it is reduced, w.r.t. the nominal atmospheric interface, is to achieve a more
     *  optimal atmosphere estimation. Below the reduced altitude the effect of bias and noise is less influential and
     *  produces higher quality results.
     */
    const double reducedAtmosphericInterfaceRadius_;

    //! Integer denoting the number of atmosphere samples required for the atmosphere estimator to be considered initiated.
    /*!
     *  Integer denoting the number of atmosphere samples required for the atmosphere estimator to be considered initiated.
     *  Once the number of atmosphere samples is reached, the atmosphere estimator uses the moving average of the last
     *  numberOfRequiredAtmosphereSamplesForInitiation_ of orbits to estimate the atmosphere model parameters.
     */
    const unsigned int numberOfRequiredAtmosphereSamplesForInitiation_;

    //! Vector denoting the direction in the body frame along which the altimeter instrument is pointing.
    const Eigen::Vector3d altimeterBodyFixedPointingDirection_;

    //! Pair denoting the lowest and highest operation altitudes of the altimeter.
    std::pair< double, double > altimeterAltitudeRange_;

    //! Pointer to the onboard environment updater.
    /*!
     *  Pointer to the onboard environment updater. The environment is updated based on the current state and time.
     *  Calling the updateEnvironment function automatically updates all dependent variables that are needed to calulate
     *  the state derivative.
     */
    boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > onboardEnvironmentUpdater_;

    //! Unordered map of states to be updated by the environment updater.
    std::unordered_map< propagators::IntegratedStateType, Eigen::VectorXd > mapOfStatesToUpdate_;

    //! Filter object to be used for estimation of state.
    boost::shared_ptr< filters::FilterBase< double, double > > navigationFilter_;

    //! Integer denoting the index of the spherical harmonics gravity exerted by the planet being orbited.
    unsigned int sphericalHarmonicsGravityIndex_;

    //! Double denoting the current time in the estimation process.
    double currentTime_;

    //! Double denoting the integration constant time step for navigation.
    double navigationRefreshStepSize_;

    //! Vector denoting the current estimated translational state in Cartesian elements.
    Eigen::Vector6d currentEstimatedCartesianState_;

    //! Vector denoting the current estimated translational state in Keplerian elements.
    Eigen::Vector6d currentEstimatedKeplerianState_;

    //! Pointer to root-finder used to esimated the time of periapsis.
    boost::shared_ptr< root_finders::BisectionCore< double > > areaBisectionRootFinder_;

    //! Pointer to integrator settings for onboard propagation.
    /*!
     *  Pointer to integrator settings for onboard propagation. Note that this object is not used to propagate the spacecraft
     *  state over time (this is carried out by the navigationFilter_ object), but is used by the navigation system or other
     *  systems to predict the spacecraft position at a future time.
     */
    boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > onboardIntegratorSettings_;

    //! Pointer to propagator settings for onboard propagation.
    /*!
     *  Pointer to propagator settings for onboard propagation. Note that this object is not used to propagate the spacecraft
     *  state over time (this is carried out by the navigationFilter_ object), but is used by the navigation system or other
     *  systems to predict the spacecraft position at a future time.
     */
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > onboardPropagatorSettings_;

    //! Enumeration denoting the previous navigation phase.
    NavigationPhaseIndicator previousNavigationPhase_;

    //! Enumeration denoting the current navigation phase.
    NavigationPhaseIndicator currentNavigationPhase_;

    //! Vector denoting the current estimated translational accelerations.
    Eigen::Vector3d currentEstimatedTranslationalAcceleration_;

    //! Vector denoting the current estimated gravitational translational accelerations.
    Eigen::Vector3d currentEstimatedGravitationalTranslationalAcceleration_;

    //! Vector denoting the current estimated non-gravitational translational accelerations.
    Eigen::Vector3d currentEstimatedNonGravitationalTranslationalAcceleration_;

    //! Vector denoting the bias and scale errors of the accelerometer after calibration.
    Eigen::Vector6d estimatedAccelerometerErrors_;

    //! History of estimated errors in Keplerian state as computed by the Periapse Time Estimator for each orbit.
    std::map< unsigned int, Eigen::Vector6d > historyOfEstimatedErrorsInKeplerianState_;

    //! Boolean denoting whether the atmosphere estimator has been initialized.
    /*!
     *  Boolean denoting whether the atmosphere estimator has been initialized. Note that the initialization is assumed to
     *  be achieved once at least numberOfRequiredAtmosphereSamplesForInitiation_ orbits have been carried out.
     */
    bool atmosphereEstimatorInitialized_;

    //! Hisotry of estimated atmosphere model parameters.
    /*!
     *  Hisotry of estimated atmosphere model parameters. These values are used as input to the moving average calculation to
     *  achieve a more accurate estimate of the atospheric properties.
     */
    std::vector< Eigen::VectorXd > historyOfEstimatedAtmosphereParameters_;

    //! History of estimated translational (in Cartesian and Keplerian elements) and rotational states.
    /*!
     *  History of estimated translational (in Cartesian and Keplerian elements) and rotational states,
     *  stored as a map, where the keys are times and the mapped values are pairs, whose first element is the pair
     *  of Cartesian and Keplerian elements, whereas the second element is the rotational state.
     */
    std::map< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > > historyOfEstimatedStates_;

    //! History of estimated translational states in Cartesian and Keplerian elements for current orbit.
    /*!
     *  History of estimated translational states in Cartesian and Keplerian elements for current orbit, stored as a map,
     *  where the keys are times and the mapped values are a pair of vectors, denoting the Cartesian and
     *  Keplerian elements.
     */
    std::map< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > > currentOrbitHistoryOfEstimatedTranslationalStates_;

    //! History of estimated translational accelerations for current orbit.
    std::map< double, Eigen::Vector3d > currentOrbitHistoryOfEstimatedTranslationalAccelerations_;

    //! History of estimated non-gravitational translational accelerations for current orbit.
    std::map< double, Eigen::Vector3d > currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_;

    //! Predefined iterator to save (de)allocation time.
    basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationMapIterator_;

    //! Predefined iterator to save (de)allocation time.
    std::map< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > >::const_iterator translationalStateIterator_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_NAVIGATION_SYSTEM_H
