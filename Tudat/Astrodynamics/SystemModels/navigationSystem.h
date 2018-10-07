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
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/aerodynamicAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicAccelerationPartial.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"
#include "Tudat/Astrodynamics/SystemModels/extraFunctions.h"
#include "Tudat/Astrodynamics/SystemModels/instrumentsModel.h"

#include "Tudat/Mathematics/Filters/createFilter.h"
#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

//! Typedefs and using statements to simplify code.
namespace Eigen { typedef Eigen::Matrix< double, 9, 1 > Vector9d; typedef Eigen::Matrix< double, 9, 9 > Matrix9d; }

namespace tudat
{

namespace system_models
{

//! Function to remove errors in inertial measurement unit measurements based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromInertialMeasurementUnitMeasurement( const Eigen::Vector3d& currentInertialMeasurementUnitMeasurement,
                                                                    const Eigen::Vector3d& inertialMeasurementUnitErrors );

//! Class for the navigation system of an aerobraking maneuver.
class NavigationSystem
{
public:

    //! Enumeration for indices of navigation state vector.
    enum NavigationStateIndices
    {
        cartesian_position_index = 0,
        cartesian_velocity_index = 3,
        accelerometer_bias_index = 6,
        quaternion_real_index = 9,
        quaternion_imaginary_index = 10,
        gyroscope_bias_index = 13,
    };

    //! Enumeration for navigation phases.
    enum NavigationPhaseIndicator
    {
        undefined_navigation_phase = -1,
        unaided_navigation_phase = 0,
        aided_navigation_phase = 1
    };

    //! Constructor.
    /*!
     *  Constructor for the navigation system of an aerobraking maneuver.
     */
    NavigationSystem( const simulation_setup::NamedBodyMap& onboardBodyMap,
                      const basic_astrodynamics::AccelerationMap& onboardAccelerationModelMap,
                      const std::string& spacecraftName, const std::string& planetName,
                      const boost::shared_ptr< filters::FilterSettings< > > navigationFilterSettings,
                      const aerodynamics::AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel,
                      const double atmosphericInterfaceAltitude, const double reducedAtmosphericInterfaceAltitude,
                      const double periapseEstimatorConstant, const unsigned int numberOfRequiredAtmosphereSamplesForInitiation,
                      const unsigned int frequencyOfDeepSpaceNetworkTracking,
                      const std::vector< Eigen::Vector3d >& altimeterPointingDirectionInAltimeterFrame = std::vector< Eigen::Vector3d >( ),
                      const reference_frames::AerodynamicsReferenceFrames altimeterFrame = reference_frames::inertial_frame,
                      const std::pair< double, double > altimeterAltitudeRange = std::make_pair( TUDAT_NAN, TUDAT_NAN ) ) :
        onboardBodyMap_( onboardBodyMap ), onboardAccelerationModelMap_( onboardAccelerationModelMap ),
        spacecraftName_( spacecraftName ), planetName_( planetName ), navigationFilterSettings_( navigationFilterSettings ),
        selectedOnboardAtmosphereModel_( selectedOnboardAtmosphereModel ),
        planetaryGravitationalParameter_( onboardBodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( ) ),
        planetaryRadius_( onboardBodyMap_.at( planetName_ )->getShapeModel( )->getAverageRadius( ) ),
        atmosphericInterfaceRadius_( planetaryRadius_ + atmosphericInterfaceAltitude ),
        reducedAtmosphericInterfaceRadius_( planetaryRadius_ + reducedAtmosphericInterfaceAltitude ),
        periapseEstimatorConstant_( periapseEstimatorConstant ),
        numberOfRequiredAtmosphereSamplesForInitiation_( numberOfRequiredAtmosphereSamplesForInitiation ),
        frequencyOfDeepSpaceNetworkTracking_( frequencyOfDeepSpaceNetworkTracking ),
        altimeterPointingDirectionInAltimeterFrame_( altimeterPointingDirectionInAltimeterFrame ),
        altimeterFrame_( altimeterFrame ), altimeterAltitudeRange_( altimeterAltitudeRange )
    {
        // Get indeces of accelerations of interest
        for ( accelerationMapConstantIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapConstantIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapConstantIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationMapConstantIterator_->second.size( ); i++ )
            {
                // Get index of spherical harmonics gravity acceleration (which is not measured by the IMUs)
                if ( ( basic_astrodynamics::getAccelerationModelType( accelerationMapConstantIterator_->second[ i ] ) ==
                       basic_astrodynamics::spherical_harmonic_gravity ) &&
                     ( accelerationMapConstantIterator_->first == planetName_ ) )
                {
                    sphericalHarmonicsGravityIndex_ = i;
                }

                // Get index of aerodynamic acceleration
                if ( ( basic_astrodynamics::getAccelerationModelType( accelerationMapConstantIterator_->second[ i ] ) ==
                       basic_astrodynamics::aerodynamic ) &&
                     ( accelerationMapConstantIterator_->first == planetName_ ) )
                {
                    aerodynamicsIndex_ = i;
                }
            }
        }

        // Create environment updater
        createOnboardEnvironmentUpdater( );

        // Set values to their initial conditions
        atmosphereEstimatorInitialized_ = false;
        estimatedAccelerometerErrors_.setZero( );
        previousNavigationPhase_ = undefined_navigation_phase;
        timeAtNavigationPhaseInterface_ = TUDAT_NAN;
    }

    //! Destructor.
    ~NavigationSystem( ) { }

    //! Function to create navigation objects for onboard state estimation.
    /*!
     *  Function to create navigation objects for onboard state estimation. This function should be called before any feature of
     *  the navigation system is used, as it creates most of the objects that are needed for state estimation (i.e., navigation
     *  filter, root-finder, onboard integrator and propagator settings).
     *  \param saveFrequency Integer denoting the frequency with which state estimates need to be stored in history.
     *  \param accelerometerMeasurementFunction Function returning the current accelerometer measurement, taken directly from the
     *      onboard inertial measurement unit.
     */
    void createNavigationSystemObjects(
            const unsigned int saveFrequency,
            const boost::function< Eigen::Vector3d( ) >& accelerometerMeasurementFunction );

    //! Function to determine the navigation phase.
    NavigationPhaseIndicator determineNavigationPhase( )
    {
        // Store previous navigation phase
        previousNavigationPhase_ = currentNavigationPhase_;

        // Declare navigation phase indicator and set value to unaided phase
        NavigationPhaseIndicator detectedNavigationPhase = unaided_navigation_phase;

        // Based on current altitude, set value to aided phase
        double currentEstimatedAltitude = currentEstimatedCartesianState_.segment( 0, 3 ).norm( ) - planetaryRadius_;
        if ( ( currentEstimatedAltitude < 7.5e6 ) || ( currentEstimatedAltitude > 25.0e6 ) )
        {
            detectedNavigationPhase = aided_navigation_phase;
        }

        // Make sure that switch only happens when it really needs to
        // Due to uncertainties in the navigation estimate, altitude may oscillate around the threshold above
        if ( previousNavigationPhase_ != detectedNavigationPhase )
        {
            // Check if navigation phase needs to be reverted
            // The reversion happens only if the phase change is detected within 5 minutes of the previous phase change
            if ( currentTime_ < ( timeAtNavigationPhaseInterface_ + ( 5.0 * 60.0 ) ) )
            {
                detectedNavigationPhase = ( detectedNavigationPhase == unaided_navigation_phase ) ? aided_navigation_phase :
                                                                                                    unaided_navigation_phase;
            }
            else
            {
                // Store time of navigation phase change
                timeAtNavigationPhaseInterface_ = currentTime_;
            }
        }

        // Set and give current navigation phase
        currentNavigationPhase_ = detectedNavigationPhase;
        return detectedNavigationPhase;
    }

    //! Function to run the State Estimator (SE).
    void runStateEstimator( const Eigen::Vector3d& currentExternalMeasurement,
                            const Eigen::Vector3d& scheduledApoapsisManeuver = Eigen::Vector3d::Zero( ) )
    {
        if ( int( ( currentTime_ - initialTime_ ) * 10.0 ) % int( 1.0e4 * 10.0 ) == 0.0 )
            std::cout << int( currentTime_ - initialTime_ ) << std::endl;

        // Add maneuver if requested
        if ( !scheduledApoapsisManeuver.isZero( ) )
        {
            Eigen::Vector6d currentEstimatedCartesianState = currentEstimatedCartesianState_;
            currentEstimatedCartesianState.segment( 3, 3 ) += scheduledApoapsisManeuver;
            setCurrentEstimatedCartesianState( currentEstimatedCartesianState );
            updateOnboardModel( );
        }

        // Perform extra functions when switching phase
        if ( ( previousNavigationPhase_ == aided_navigation_phase ) && ( currentNavigationPhase_ == unaided_navigation_phase ) )
        {
            // Improve state estimate if passing from aided to unaided
            improveStateEstimateOnNavigationPhaseTransition( );
        }
        else if ( ( previousNavigationPhase_ == unaided_navigation_phase ) && ( currentNavigationPhase_ == aided_navigation_phase ) )
        {
            // Reset covariance matrix if passing from unaided to aided
            Eigen::Vector9d currentNavigationFilterState = navigationFilter_->getCurrentStateEstimate( );
            Eigen::Matrix9d currentNavigationFilterCovarianceMatrix = navigationFilter_->getCurrentCovarianceEstimate( );
            currentNavigationFilterCovarianceMatrix.setIdentity( );
            navigationFilter_->modifyCurrentStateAndCovarianceEstimates( currentNavigationFilterState,
                                                                         currentNavigationFilterCovarianceMatrix );
        }

        // Update current state based on the detected navigation phase
        switch ( currentNavigationPhase_ )
        {
        case unaided_navigation_phase:
        {
            // Update Cartesian state
            Eigen::Vector6d currentEstimatedCartesianState = currentEstimatedCartesianState_;
            Eigen::Vector6d currentEstimatedCartesianStateDerivative;
            currentEstimatedCartesianStateDerivative.segment( 0, 3 ) = currentEstimatedCartesianState.segment( 3, 3 );
            currentEstimatedCartesianStateDerivative.segment( 3, 3 ) = currentEstimatedTranslationalAcceleration_;
            currentEstimatedCartesianState += navigationRefreshStepSize_ * currentEstimatedCartesianStateDerivative;

            // Update time
            currentTime_ += navigationRefreshStepSize_;

            // Update navigation system
            setCurrentEstimatedCartesianState( currentEstimatedCartesianState );
            break;
        }
        case aided_navigation_phase:
        {
            // Update filter
            navigationFilter_->updateFilter( currentExternalMeasurement );

            // Extract tiem and estimated state and update navigation system
            currentTime_ = navigationFilter_->getCurrentTime( );
            Eigen::Vector9d updatedEstimatedState = navigationFilter_->getCurrentStateEstimate( );
            setCurrentEstimatedCartesianState( updatedEstimatedState.segment( 0, 6 ) );
            break;
        }
        default:
            throw std::runtime_error( "Error in navigation system. Current navigation (" +
                                      std::to_string( currentNavigationPhase_ ) + ") phase not supported." );
        }

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
        for ( translationalStateConstantIterator_ = currentOrbitHistoryOfEstimatedTranslationalStates_.begin( );
              translationalStateConstantIterator_ != currentOrbitHistoryOfEstimatedTranslationalStates_.end( );
              translationalStateConstantIterator_++ )
        {
            if ( translationalStateConstantIterator_->second.first.segment( 0, 3 ).norm( ) <= atmosphericInterfaceRadius_ )
            {
                // Retireve time, state and acceleration of where the altitude is below the atmospheric interface
                currentIterationTime = translationalStateConstantIterator_->first;
                mapOfEstimatedCartesianStatesBelowAtmosphericInterface[ currentIterationTime ] =
                        translationalStateConstantIterator_->second.first;
                mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ] =
                        translationalStateConstantIterator_->second.second;
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
            postProcessAccelerometerMeasurements( vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface );

            // Compute magnitude of aerodynamic acceleration
            std::vector< double > vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface;
            for ( unsigned int i = 0; i < vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
            {
                vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back(
                            vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ).norm( ) );
            }

            // Run periapse time estimator if ... (TBD)
            if ( historyOfEstimatedAtmosphereParameters_.size( ) > 0 )
            {
                runPeriapseTimeEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                          vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
            }

            // Run atmosphere estimator with processed results
            runAtmosphereEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
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
        // Inform user
        std::cout << std::endl << "Processing Deep Space Network Tracking Information." << std::endl;

        // Split Deep Space Network tracking data in light-time delay and Cartesian state
        double currentLightTimeDelay = currentDeepSpaceNetworkTrackingData.first;
        Eigen::Vector6d currentTrackedState = currentDeepSpaceNetworkTrackingData.second;

        // Propagate tracked state to current time
        boost::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings =
                boost::make_shared< propagators::PropagationTimeTerminationSettings >( currentTime_ );
        Eigen::Vector6d propagatedStateBasedOnTracking = propagateTranslationalStateWithCustomTerminationSettings(
                    terminationSettings, currentTrackedState, currentTime_ - currentLightTimeDelay ).second.first.rbegin( )->second;

        // Convert to Keplerian elements and reset true anomaly
        Eigen::Vector6d propagatedStateBasedOnTrackingInKeplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                    propagatedStateBasedOnTracking, planetaryGravitationalParameter_ );
        propagatedStateBasedOnTrackingInKeplerianElements[ 5 ] = currentEstimatedKeplerianState_[ 5 ]; // reset true anomaly

        // Inform user
        std::cout << "Propagated state for " << currentLightTimeDelay / 60.0 << " minutes." << std::endl;

        // Reset navigation filter (including covariance)
        Eigen::Matrix9d currentEstimatedCovarianceMatrix = navigationFilter_->getCurrentCovarianceEstimate( );
        currentEstimatedCovarianceMatrix.block( 0, 0, 6, 6 ).setIdentity( );
        setCurrentEstimatedKeplerianState( propagatedStateBasedOnTrackingInKeplerianElements, currentEstimatedCovarianceMatrix );
    }

    //! Function to propagate translational Cartesian state to specified termination settings.
    /*!
     *  Function to propagate translational Cartesian state to specified termination settings. Note that this function does not
     *  continuously estimate the spacecraft state, but instead it is called by the navigation system or any other onboard system
     *  to propagate the current estimated state with the onboard model to some termination condition. The initial time and state are
     *  automatically set to the current time and current estimated translational Cartesian state, respectively, but they can be
     *  overwritten by setting them as inputs.
     *  \param propagationTerminationSettings Termination settings to be used to stop the propagation.
     *  \param initialTranslationalCartesianState Vector denoting the initial conditions for the translational state in Cartesian elements.
     *  \param initialTime Double denoting the initial time for propagation.
     *  \return Pair of boolean denoting whether propagation was successful and pair of dynamics simulator results. The first element
     *      of the second pair is the map of time and propagated state, whereas the second one is the map of time and dependent
     *      variable values.
     */
    std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > >
    propagateTranslationalStateWithCustomTerminationSettings(
            const boost::shared_ptr< propagators::PropagationTerminationSettings > propagationTerminationSettings,
            const Eigen::Vector6d& initialTranslationalCartesianState = Eigen::Vector6d::Zero( ),
            const double initialTime = -1.0 )
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
        if ( initialTranslationalCartesianState.isZero( ) )
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
        bool wasPropagationSuccessful = true;
        propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
                    onboardBodyMap_, onboardIntegratorSettings_, onboardPropagatorSettings_ );
        if ( dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) !=
             propagators::termination_condition_reached )
        {
            wasPropagationSuccessful = false;
        }

        // Retrieve results from onboard computer and systems
        std::map< double, Eigen::VectorXd > translationalStateResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariablesResult = dynamicsSimulator.getDependentVariableHistory( );

        // Give output
        return std::make_pair( wasPropagationSuccessful, std::make_pair( translationalStateResult, dependentVariablesResult ) );
    }

    //! Function to retireve current time.
    double getCurrentTime( ) { return currentTime_; }

    //! Function to retireve current state.
    Eigen::Vector9d getCurrentEstimatedState( ) { return navigationFilter_->getCurrentStateEstimate( ); }

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
    double getDensityAtSpecifiedConditions( const Eigen::Vector6d& estimatedState = Eigen::Vector6d::Zero( ) )
    {
        if ( estimatedState.isZero( ) )
        {
            return onboardBodyMap_.at( planetName_ )->getAtmosphereModel( )->getDensity(
                        onboardBodyMap_.at( spacecraftName_ )->getFlightConditions( )->getCurrentAltitude( ),
                        onboardBodyMap_.at( spacecraftName_ )->getFlightConditions( )->getCurrentLongitude( ), 0.0, 0.0 );
        }
        else
        {
            return onboardBodyMap_.at( planetName_ )->getAtmosphereModel( )->getDensity(
                        estimatedState.segment( 0, 3 ).norm( ) - planetaryRadius_, 0.0, 0.0, 0.0 );
        }
    }

    //! Function to retrieve whether the spacecraft altitude is above the dynamic atmospheric interface line.
    bool getIsSpacecraftAboveDynamicAtmosphericInterfaceAltitude( )
    {
        // Compute current DAIA
        double currentDynamicAtmosphericInterfaceRadius = 0.25 * ( currentEstimatedKeplerianState_[ 0 ] *
                ( 1.0 + currentEstimatedKeplerianState_[ 1 ] ) ) + planetaryRadius_;

        // Compute current radial distance
        double currentRadialDistance = currentEstimatedCartesianState_.segment( 0, 3 ).norm( );

        // Output whether the altitude is below DAIA
        return ( currentRadialDistance > currentDynamicAtmosphericInterfaceRadius );
    }

    //! Function to retrieve current estimated translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated translational accelerations exerted on the spacecraft. The acceleration
     *  is computed by using the onboard accelerations model map.
     *  \param estimatedState State at which the environment has to be updated and accelerations computed and stored. If left
     *      empty, the accelerations at currentEstimatedCartesianState_ are retrieved.
     *  \return Vector denoting the full acceleration vector experienced by the spacecraft.
     */
    Eigen::Vector3d getCurrentEstimatedTranslationalAcceleration(
            const Eigen::Vector6d& estimatedState = Eigen::Vector6d::Zero( ) )
    {
        if ( !estimatedState.isZero( ) )
        {
            if ( !estimatedState.isApprox( currentEstimatedCartesianState_ ) )
            {
                updateOnboardModel( estimatedState );
            }
        }
        return currentEstimatedTranslationalAcceleration_;
    }

    //! Function to retrieve current estimated gravitational translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated gravitational translational accelerations exerted on the spacecraft. The
     *  acceleration is computed by using the onboard accelerations model map. Note that this vector represents only the
     *  gravitational accelerations.
     *  \param estimatedState State at which the environment has to be updated and accelerations computed and stored. If left
     *      empty, the accelerations at currentEstimatedCartesianState_ are retrieved.
     *  \return Vector denoting only the gravitational accelerations experienced by the spacecraft.
     */
    Eigen::Vector3d getCurrentEstimatedGravitationalTranslationalAcceleration(
            const Eigen::Vector6d& estimatedState = Eigen::Vector6d::Zero( ) )
    {
        if ( !estimatedState.isZero( ) )
        {
            if ( !estimatedState.isApprox( currentEstimatedCartesianState_ ) )
            {
                updateOnboardModel( estimatedState );
            }
        }
        return currentEstimatedGravitationalTranslationalAcceleration_;
    }

    //! Function to retrieve current estimated non-gravitational translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated non-gravitational translational accelerations exerted on the spacecraft. The
     *  acceleration is computed by using the onboard accelerations model map. Note that this vector represents only the
     *  non-gravitational accelerations (which are supposed to emulate the accelerations measured by the IMU).
     *  \param estimatedState State at which the environment has to be updated and accelerations computed and stored. If left
     *      empty, the accelerations at currentEstimatedCartesianState_ are retrieved.
     *  \return Vector denoting only the non-gravitational (i.e., aerodynamic) accelerations experienced by the spacecraft.
     */
    Eigen::Vector3d getCurrentEstimatedNonGravitationalTranslationalAcceleration(
            const Eigen::Vector6d& estimatedState = Eigen::Vector6d::Zero( ) )
    {
        if ( !estimatedState.isZero( ) )
        {
            if ( !estimatedState.isApprox( currentEstimatedCartesianState_ ) )
            {
                updateOnboardModel( estimatedState );
            }
        }
        return currentEstimatedNonGravitationalTranslationalAcceleration_;
    }

    //! Function to retrieve the estimated accelerometer errors.
    Eigen::Vector3d getEstimatedAccelerometerErrors( )
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

    //! Function to retrieve history of estimated atmosphere parameters.
    std::map< unsigned int, Eigen::VectorXd > getHistoryOfEstimatedAtmosphereParameters( )
    {
        return historyOfEstimatedAtmosphereParameters_;
    }

    //! Function to retrieve history of estimated translational accelerations for the current orbit.
    std::map< double, Eigen::Vector3d > getCurrentOrbitHistoryOfEstimatedTranslationalAccelerations( )
    {
        return currentOrbitHistoryOfEstimatedTranslationalAccelerations_;
    }

    //! Function to retrieve history of estimated gravitational translational accelerations for the current orbit.
    std::map< double, Eigen::Vector3d > getCurrentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations( )
    {
        return currentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations_;
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

    //! Function to retireve the frequency (in days) with which Deep Space Network tracking has to be performed.
    unsigned int getFrequencyOfDeepSpaceNetworkTracking( ) { return frequencyOfDeepSpaceNetworkTracking_; }

    void setCurrentTime( const double newCurrentTime )
    {
        currentTime_ = newCurrentTime;
        navigationFilter_->modifyCurrentTime( newCurrentTime );
    }

    //! Function to set current Cartesian state to new value.
    /*!
     *  Function to set current Cartesian state to new value. The value of the Keplerian state is set
     *  automatically by converting the Cartesian state to Keplerian elements. Also, the function updates the
     *  history of translational to the new values.
     *  \param newCurrentKeplerianState New Cartesian state at current time.
     *  \param newCurrentCovarianceMatrix New covariance matrix at current time.
     */
    void setCurrentEstimatedCartesianState( const Eigen::Vector6d& newCurrentCartesianState,
                                            const Eigen::Matrix9d& newCurrentCovarianceMatrix = Eigen::Matrix9d::Zero( ) )
    {
        // Update navigation system current value
        currentEstimatedCartesianState_ = newCurrentCartesianState;
        currentEstimatedKeplerianState_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                    newCurrentCartesianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values

        // Update navigation filter current value
        Eigen::Vector9d updatedCurrentEstimatedState = navigationFilter_->getCurrentStateEstimate( );
        updatedCurrentEstimatedState.segment( 0, 6 ) = currentEstimatedCartesianState_;
        navigationFilter_->modifyCurrentStateAndCovarianceEstimates( updatedCurrentEstimatedState,
                                                                     newCurrentCovarianceMatrix );

        // Update time in filter if Kepler phase
        if ( currentNavigationPhase_ == unaided_navigation_phase )
        {
            navigationFilter_->modifyCurrentTime( currentTime_ );
        }
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
                                            const Eigen::Matrix9d& newCurrentCovarianceMatrix = Eigen::Matrix9d::Zero( ) )
    {
        // Update navigation system current value
        currentEstimatedKeplerianState_ = newCurrentKeplerianState;
        currentEstimatedCartesianState_ = orbital_element_conversions::convertKeplerianToCartesianElements(
                    newCurrentKeplerianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values

        // Update navigation filter current value
        Eigen::Vector9d updatedCurrentEstimatedState = navigationFilter_->getCurrentStateEstimate( );
        updatedCurrentEstimatedState.segment( 0, 6 ) = currentEstimatedCartesianState_;
        navigationFilter_->modifyCurrentStateAndCovarianceEstimates( updatedCurrentEstimatedState,
                                                                     newCurrentCovarianceMatrix );

        // Update time in filter if Kepler phase
        if ( currentNavigationPhase_ == unaided_navigation_phase )
        {
            navigationFilter_->modifyCurrentTime( currentTime_ );
        }
    }

    //! Function to set the apoapsis value of the Keplerian state.
    void setEstimatedApoapsisKeplerianState( )
    {
        estimatedKeplerianStateAtPreviousApoapsis_ = currentEstimatedKeplerianState_;
    }

    //! Reset navigation filter integration step size.
    void resetNavigationRefreshStepSize( const double newNavigationRefreshStepSize )
    {
        atmosphericNavigationRefreshStepSize_ = std::min( navigationRefreshStepSize_, newNavigationRefreshStepSize );
        navigationRefreshStepSize_ = newNavigationRefreshStepSize;
        navigationFilter_->modifyFilteringStepSize( newNavigationRefreshStepSize );
    }

    //! Function to clear history of estimated states and accelerations for the current orbit.
    void clearCurrentOrbitEstimationHistory( )
    {
        // Empty current orbit history maps
        currentOrbitHistoryOfEstimatedTranslationalStates_.clear( );
        currentOrbitHistoryOfEstimatedTranslationalAccelerations_.clear( );
        currentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations_.clear( );
        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_.clear( );

        // Empty filter estimates
        navigationFilter_->clearFilterHistory( );
    }

    //! Function to revert to the previous time step.
    /*!
     *  Function to revert to the previous time step. This function is run if the current propagation needs to be stopped, since
     *  the current time will be run the next time the GNC system is called.
     *  \param currentTime Double denoting the current time, i.e., the instant that has to be discarded.
     */
    void revertToPreviousTimeStep( const double currentTime )
    {
        // Erase the entry corresponding to the input time for the current orbit history maps
        if ( currentOrbitHistoryOfEstimatedTranslationalStates_.count( currentTime ) != 0 )
        {
            currentOrbitHistoryOfEstimatedTranslationalStates_.erase( currentTime );
        }
        if ( currentOrbitHistoryOfEstimatedTranslationalAccelerations_.count( currentTime ) != 0 )
        {
            currentOrbitHistoryOfEstimatedTranslationalAccelerations_.erase( currentTime );
        }
        if ( currentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations_.count( currentTime ) != 0 )
        {
            currentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations_.erase( currentTime );
        }
        if ( currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_.count( currentTime ) != 0 )
        {
            currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_.erase( currentTime );
        }

        // Empty filter estimates
        navigationFilter_->revertToPreviousTimeStep( currentTime );

        // Revert to previous estimates
        currentTime_ = currentOrbitHistoryOfEstimatedTranslationalStates_.rbegin( )->first;
        currentEstimatedCartesianState_ = currentOrbitHistoryOfEstimatedTranslationalStates_.rbegin( )->second.first;
        currentEstimatedKeplerianState_ = currentOrbitHistoryOfEstimatedTranslationalStates_.rbegin( )->second.second;
        currentEstimatedTranslationalAcceleration_ = currentOrbitHistoryOfEstimatedTranslationalAccelerations_.rbegin( )->second;
        currentEstimatedGravitationalTranslationalAcceleration_ =
                currentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations_.rbegin( )->second;
        currentEstimatedNonGravitationalTranslationalAcceleration_ =
                currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_.rbegin( )->second;
    }

    //! Integer denoting the current orbit counter.
    unsigned int currentOrbitCounter_;

private:

    //! Function to create the onboard environment updater.
    void createOnboardEnvironmentUpdater( );

    //! Function to update the body and acceleration map with the current time and state information.
    /*!
     *  Function to update the body and acceleration map with the current time and state information.
     *  \param estimatedState State at which the environment has to be updated and accelerations computed and stored. If left
     *      empty, currentEstimatedCartesianState_ is taken as state.
     */
    void updateOnboardModel( const Eigen::Vector6d& estimatedState = Eigen::Vector6d::Zero( ) )
    {
        // Update environment
        if ( estimatedState.isZero( ) )
        {
            mapOfStatesToUpdate_[ propagators::translational_state ] = currentEstimatedCartesianState_;
        }
        else
        {
            mapOfStatesToUpdate_[ propagators::translational_state ] = estimatedState;
        }
        mapOfStatesToUpdate_[ propagators::translational_state ] += onboardBodyMap_.at( planetName_ )->getState( );
        onboardEnvironmentUpdater_->updateEnvironment( currentTime_, mapOfStatesToUpdate_ );

        // Loop over bodies exerting accelerations on spacecraft
        for ( accelerationMapConstantIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapConstantIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapConstantIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationMapConstantIterator_->second.size( ); i++ )
            {
                // Update acceleration model
                accelerationMapConstantIterator_->second[ i ]->updateMembers( currentTime_ );
            }
        }

        // Update Jacobians based on acceleration partials
        onboardGravitationalAccelerationPartials_->update( currentTime_ );
        onboardAerodynamicAccelerationPartials_->update( currentTime_ );

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
        for ( accelerationMapConstantIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapConstantIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapConstantIterator_++ )
        {
            // Loop over each accelerations
            for ( unsigned int i = 0; i < accelerationMapConstantIterator_->second.size( ); i++ )
            {
                // Calculate acceleration and add to state derivative
                currentAcceleration = accelerationMapConstantIterator_->second[ i ]->getAcceleration( );
                currentEstimatedTranslationalAcceleration_ += currentAcceleration;

                // Only add the gravitational accelerations
                if ( ( i == sphericalHarmonicsGravityIndex_ ) && ( accelerationMapConstantIterator_->first == planetName_ ) )
                {
                    currentEstimatedGravitationalTranslationalAcceleration_ += currentAcceleration;
                }
                // Disregard the gravitational acceleration, since IMUs do not measure them
                else
                {
                    currentEstimatedNonGravitationalTranslationalAcceleration_ += currentAcceleration;
                }
            }
        }

        // Store acceleration value
        currentOrbitHistoryOfEstimatedTranslationalAccelerations_[ currentTime_ ] = currentEstimatedTranslationalAcceleration_;
        currentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations_[ currentTime_ ] =
                currentEstimatedGravitationalTranslationalAcceleration_;
        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_[ currentTime_ ] =
                currentEstimatedNonGravitationalTranslationalAcceleration_;
    }

    //! Function to improve the estimate of the translational state when transitioning from aided to unaided navigation.
    /*!
     *  Function to improve the estimate of the translational state when transitioning from aided to unaided navigation.
     *  The improvement algorithm uses the Keplerian elements, since they show the least variation with time. The first five
     *  elements are improved by taking the median of the last 60 seconds of estimation, whereas the true anomaly is improved by
     *  fitting a simple parabolic function through the true anomaly history, and computing the value for the current time.
     */
    void improveStateEstimateOnNavigationPhaseTransition( );

    //! Function to post-process the accelerometer measurements.
    /*!
     *  Function to post-process the accelerometer measurements. This function first applies a simple smoothing method
     *  (moving average) to the filter, then estimates the errors in the accelerometer (bias and scale factor) and then
     *  transforms the results from body-fixed frame (the frame of the accelerometer) to inertial (the frame in which
     *  accelerations are needed for propagation).
     *  \param vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface Vector of estimated aerodynamic acceleration
     *      below the atmospheric interface altitude and corrupted by accelerometer errors. The acceleration is thus taken
     *      directly from the IMU.
     */
    void postProcessAccelerometerMeasurements(
            std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface );

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
            std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
            const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

    //! Function to run the Atmosphere Estimator (AE).
    /*!
     *  Function to estimate the atmospheric properties based on the selected onboard atmosphere model.
     *  \param mapOfEstimatedKeplerianStatesBelowAtmosphericInterface Map of time and estimated Keplerian elements below
     *      the atmospheric interface altitude.
     *  \param mapOfMeasuredAerodynamicAcceleration Map of time and estimated aerodynamic acceleration. The acceleration is
     *      the one computed from the IMU measurements and is stored as a three-dimensional vector.
     */
    void runAtmosphereEstimator(
            const std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
            const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

    //! Function to model the onboard system dynamics based on the simplified onboard model.
    Eigen::Vector9d onboardSystemModel( const double currentTime, const Eigen::Vector9d& currentEstimatedState,
                                        const boost::function< Eigen::Vector3d( ) >& accelerometerMeasurementFunction );

    //! Function to model the onboard measurements based on the simplified onboard model.
    Eigen::Vector3d onboardMeasurementModel( const double currentTime, const Eigen::Vector9d& currentEstimatedState );

    //! Function to model the onboard system Jacobian based on the simplified onboard model.
    Eigen::Matrix9d onboardSystemJacobian( const double currentTime, const Eigen::Vector9d& currentEstimatedState );

    //! Function to model the onboard measurements Jacobian based on the simplified onboard model.
    Eigen::Matrix< double, 3, 9 > onboardMeasurementJacobian( const double currentTime, const Eigen::Vector9d& currentEstimatedState );

    //! Function to store current time and current state estimates.
    void storeCurrentTimeAndStateEstimates( )
    {
        // Store translational values for current orbit
        currentOrbitHistoryOfEstimatedTranslationalStates_[ currentTime_ ] = std::make_pair( currentEstimatedCartesianState_,
                                                                                             currentEstimatedKeplerianState_ );

        // Save estimated state to full history (depending on save index)
        saveIndex_ = saveIndex_ % saveFrequency_;
        if ( saveIndex_ == 0 )
        {
            historyOfEstimatedStates_[ currentTime_ ] = currentOrbitHistoryOfEstimatedTranslationalStates_[ currentTime_ ];
        }
        saveIndex_++;
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

    //! Double denoting the multiplier to account for non-Keplerian orbit.
    /*!
     *  Double denoting the multiplier to account for non-Keplerian orbit. It is used to multiply the value of estimated change
     *  in velocity of the Periapse Time Estimator, to remove the error of the assumption of Kepler orbit and impulsive drag.
     */
    const double periapseEstimatorConstant_;

    //! Integer denoting the number of atmosphere samples required for the atmosphere estimator to be considered initiated.
    /*!
     *  Integer denoting the number of atmosphere samples required for the atmosphere estimator to be considered initiated.
     *  Once the number of atmosphere samples is reached, the atmosphere estimator uses the moving average of the last
     *  numberOfRequiredAtmosphereSamplesForInitiation_ of orbits to estimate the atmosphere model parameters.
     */
    const unsigned int numberOfRequiredAtmosphereSamplesForInitiation_;

    //! Integer denoting the number of frequency of Deep Space Network tracking in days.
    /*!
     *  Integer denoting the number of frequency of Deep Space Network tracking in days. This means that every
     *  frequencyOfDeepSpaceNetworkTracking_ days the Deep Space Network tracking information will be processed and used to update
     *  the current estimated state.
     */
    const unsigned int frequencyOfDeepSpaceNetworkTracking_;

    //! Vector denoting the directions in the altimeter frame along which the altimeter instruments are pointing.
    const std::vector< Eigen::Vector3d > altimeterPointingDirectionInAltimeterFrame_;

    //! Enumeration denoting the frame in which the altimeter pointing directions are defined;
    const reference_frames::AerodynamicsReferenceFrames altimeterFrame_;

    //! Pair denoting the lowest and highest operation altitudes of the altimeter.
    const std::pair< double, double > altimeterAltitudeRange_;

    //! Integer denoting the frequency with which state estimates need to be stored in history.
    unsigned int saveFrequency_;

    //! Integer denoting the current save index.
    unsigned int saveIndex_;

    //! Unordered map of states to be updated by the environment updater.
    std::unordered_map< propagators::IntegratedStateType, Eigen::VectorXd > mapOfStatesToUpdate_;

    //! Pointer to the onboard environment updater.
    /*!
     *  Pointer to the onboard environment updater. The environment is updated based on the current state and time.
     *  Calling the updateEnvironment function automatically updates all dependent variables that are needed to calulate
     *  the state derivative.
     */
    boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > onboardEnvironmentUpdater_;

    //! Pointer to the onboard spherical harmonics gravity acceleration partial object.
    /*!
     *  Pointer to the onboard spherical harmonics gravity acceleration partial object, which is used to compute the Jacobian of
     *  the spherical harmonics acceleration w.r.t. the current Cartesian state.
     */
    boost::shared_ptr< acceleration_partials::SphericalHarmonicsGravityPartial > onboardGravitationalAccelerationPartials_;

    //! Pointer to the onboard aerodynamic acceleration partial object.
    /*!
     *  Pointer to the onboard aerodynamic acceleration partial object, which is used to compute the Jacobian of
     *  the aerodynamics acceleration w.r.t. the current Cartesian state.
     */
    boost::shared_ptr< acceleration_partials::AerodynamicAccelerationPartial > onboardAerodynamicAccelerationPartials_;

    //! Filter object to be used for estimation of state.
    boost::shared_ptr< filters::FilterBase< double, double > > navigationFilter_;

    //! Integer denoting the index of the spherical harmonics gravity acceleration exerted by the planet being orbited.
    unsigned int sphericalHarmonicsGravityIndex_;

    //! Integer denoting the index of the aerodynamic acceleration exerted by the planet being orbited.
    unsigned int aerodynamicsIndex_;

    //! Double denoting the initial time in the estimation process.
    double initialTime_;

    //! Double denoting the current time in the estimation process.
    double currentTime_;

    //! Double denoting the integration constant time step for navigation.
    double navigationRefreshStepSize_;

    //! Double denoting the integration constant time step for navigation during the atmospheric phase.
    double atmosphericNavigationRefreshStepSize_;

    //! Vector denoting the current estimated translational state in Cartesian elements.
    Eigen::Vector6d currentEstimatedCartesianState_;

    //! Vector denoting the current estimated translational state in Keplerian elements.
    Eigen::Vector6d currentEstimatedKeplerianState_;

    //! Vector denoting Keplerian elements at the previous apoapsis.
    /*!
     *  Vector denoting Keplerian elements at the previous apoapsis. This vector is used by the Periapse Time Estimator to extract the
     *  information on semi-major axis and eccentricity.
     */
    Eigen::Vector6d estimatedKeplerianStateAtPreviousApoapsis_;

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

    //! Double denoting the time at which the navigation phase has changed.
    double timeAtNavigationPhaseInterface_;

    //! Vector denoting the current estimated translational accelerations.
    Eigen::Vector3d currentEstimatedTranslationalAcceleration_;

    //! Vector denoting the current estimated gravitational translational accelerations.
    Eigen::Vector3d currentEstimatedGravitationalTranslationalAcceleration_;

    //! Vector denoting the current estimated non-gravitational translational accelerations.
    Eigen::Vector3d currentEstimatedNonGravitationalTranslationalAcceleration_;

    //! Vector denoting the bias and scale errors of the accelerometer after calibration.
    Eigen::Vector3d estimatedAccelerometerErrors_;

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
    std::map< unsigned int, Eigen::VectorXd > historyOfEstimatedAtmosphereParameters_;

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

    //! History of estimated gravitational translational accelerations for current orbit.
    std::map< double, Eigen::Vector3d > currentOrbitHistoryOfEstimatedGravitationalTranslationalAccelerations_;

    //! History of estimated non-gravitational translational accelerations for current orbit.
    std::map< double, Eigen::Vector3d > currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_;

    //! Predefined iterator to save (de)allocation time.
    basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationMapConstantIterator_;

    //! Predefined iterator to save (de)allocation time.
    std::map< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > >::const_iterator translationalStateConstantIterator_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_NAVIGATION_SYSTEM_H
