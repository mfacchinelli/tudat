/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Facchinelli, M. (2018). Aerobraking Navigation, Guidance and Control.
 *          Master Thesis, Delft University of Technology.
 */

#ifndef TUDAT_ONBOARD_COMPUTER_MODEL_H
#define TUDAT_ONBOARD_COMPUTER_MODEL_H

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/Propagators/imanRmsMethod.h"

#include "Tudat/Astrodynamics/SystemModels/controlSystem.h"
#include "Tudat/Astrodynamics/SystemModels/guidanceSystem.h"
#include "Tudat/Astrodynamics/SystemModels/instrumentsModel.h"
#include "Tudat/Astrodynamics/SystemModels/navigationSystem.h"

//! Typedefs and using statements to simplify code.
namespace Eigen { typedef Eigen::Matrix< double, 12, 1 > Vector12d; }

namespace tudat
{

namespace system_models
{

//! Class for the onboard computer of the spacecraft.
class OnboardComputerModel
{
public:

    //! Constructor.
    OnboardComputerModel( const boost::shared_ptr< ControlSystem > controlSystem,
                          const boost::shared_ptr< GuidanceSystem > guidanceSystem,
                          const boost::shared_ptr< NavigationSystem > navigationSystem,
                          const boost::shared_ptr< InstrumentsModel > instrumentsModel,
                          const unsigned int saveFrequency = 1 ) :
        controlSystem_( controlSystem ), guidanceSystem_( guidanceSystem ), navigationSystem_( navigationSystem ),
        instrumentsModel_( instrumentsModel )
    {
        // Define internal variables
        maneuveringPhaseComplete_ = false; // simulation starts at apoapsis with possible need to perform a maneuver
        atmosphericPhaseComplete_ = true;

        performManeuverOnNextCall_ = false;
        deepSpaceNetworkTrackingInformation_ = std::make_pair( false, static_cast< unsigned int >( -1 ) );

        // Create navigation system objects
        navigationSystem_->createNavigationSystemObjects( saveFrequency );
        initialTime_ = navigationSystem_->getCurrentTime( );

        // Create guidance system objects
        guidanceSystem_->setCurrentOrbitCounter( navigationSystem_->currentOrbitCounter_ );
        guidanceSystem_->createGuidanceSystemObjects(
                    boost::bind( &NavigationSystem::propagateTranslationalStateWithCustomTerminationSettings,
                                 navigationSystem_, _1, _2, -1.0 ),
                    navigationSystem_->getStandardGravitationalParameter( ), navigationSystem_->getRadius( ) );
    }

    //! Destructor.
    ~OnboardComputerModel( ) { }

    //! Function to check whether the propagation is to be be stopped.
    /*!
     *  Function to check whether the propagation is to be be stopped, based on the estimated state. The propagation
     *  is stopped if the pacecraft is estimated to be at apoapsis. In case the spacecraft is estimated to have just
     *  left the atmosphere, the post-atmosphere processes are initiated. To make sure that the propagation is checked
     *  every time step, this function needs to be fed to the propagator via the custom propagation termination settings.
     *  \param currentTime Current time in propagation.
     *  \return Boolean denoting whether the propagation has to be stopped.
     */
    bool checkStopCondition( const double currentTime )
    {
        using mathematical_constants::PI;

        // Define output value
        bool isPropagationToBeStopped = false;

        // Check whether current step has already been performed
        if ( ( currentTime - ( navigationSystem_->getCurrentTime( ) + navigationSystem_->getNavigationRefreshStepSize( ) ) ) <= 1e-5 )
        {
            // Update instrument models
            instrumentsModel_->updateInstruments( currentTime );

            // Extract measurements
            Eigen::Vector3d currentExternalMeasurement = instrumentsModel_->getCurrentAccelerometerMeasurement( );

            // Update filter to current time
            NavigationSystem::NavigationPhaseIndicator currentNavigationPhase = navigationSystem_->determineNavigationPhase( );
            if ( performManeuverOnNextCall_ )
            {
                // Feed maneuver to the navigation system and update filter
                performManeuverOnNextCall_ = false; // reset flag
                navigationSystem_->runStateEstimator( currentExternalMeasurement,
                                                      controlSystem_->getScheduledApsisManeuver( ) );
            }
            else
            {
                // Update filter only
                navigationSystem_->runStateEstimator( currentExternalMeasurement );
            }

            // Check if it is time for a Deep Space Network update
            // The Deep Space Network tracking is scheduled every N days (where N comes from the function getFrequencyOfDeepSpaceNetworkTracking
            // of the navigation system). The following if statement makes sure that the processing of the measurements is carried out at this
            // frequency and that it is not carried out on the very first apoapsis.
            unsigned int currentDay = static_cast< unsigned int >( ( currentTime - initialTime_ ) / physical_constants::JULIAN_DAY );
            if ( ( ( currentDay % navigationSystem_->getFrequencyOfDeepSpaceNetworkTracking( ) ) == 0 ) &&
                 !deepSpaceNetworkTrackingInformation_.first && ( currentDay != deepSpaceNetworkTrackingInformation_.second ) &&
                 ( navigationSystem_->currentOrbitCounter_ > 0 ) )
            {
                deepSpaceNetworkTrackingInformation_.first = true;
                deepSpaceNetworkTrackingInformation_.second = currentDay;
            }

            // Access current state and extract true anomaly
            std::pair< Eigen::Vector6d, Eigen::Vector6d > currentEstimatedState = navigationSystem_->getCurrentEstimatedTranslationalState( );
            double currentEstimatedTrueAnomaly = currentEstimatedState.second[ 5 ];

            // Check which orbit phase the spacecraft is in (if any)
            // Check true anomaly to see if apoapsis maneuvering phase
            if ( ( currentEstimatedTrueAnomaly >= PI ) && !maneuveringPhaseComplete_ )
            {
                // Inform user
                std::cout << std::endl << "REACHED APOAPSIS. Preparing to perform maneuver." << std::endl;

                // Process Deep Space Network tracking data
                if ( deepSpaceNetworkTrackingInformation_.first )
                {
                    deepSpaceNetworkTrackingInformation_.first = false;
                    navigationSystem_->processDeepSpaceNetworkTracking( instrumentsModel_->getCurrentDeepSpaceNetworkMeasurement( ) );
                    currentEstimatedState = navigationSystem_->getCurrentEstimatedTranslationalState( ); // overwrite current state
                }

                // Store new value of apoapsis Keplerian state
                navigationSystem_->setEstimatedApoapsisKeplerianState( );

                // Determine in which phase of aerobraking the spacecraft is and perform corridor estimation
                guidanceSystem_->determineAerobrakingPhase( currentEstimatedState.second,
                                                            navigationSystem_->getAtmosphereInitiationIndicators( ) );
                guidanceSystem_->runCorridorEstimator( currentTime, currentEstimatedState.first, currentEstimatedState.second );

                // Check whether propagation needs to be stopped
                // Propagation is stopped only if an apoapsis maneuver needs to be performed
                isPropagationToBeStopped = guidanceSystem_->getIsApoapsisManeuverToBePerformed( );

                if ( propagators::IMAN_ANALYSIS_INDEX == 2 )
                {
                    std::cerr << "Apoapsis maneuver is OFF." << std::endl;
                    isPropagationToBeStopped = false;
                }

                // If propagation is to be stopped, to add Delta V to state, run remaining pre-maneuver processes
                if ( isPropagationToBeStopped )
                {
                    // Set flag to add maneuver on next function call
                    performManeuverOnNextCall_ = true;

                    // Run maneuver estimator
                    guidanceSystem_->runApoapsisManeuverEstimator( currentEstimatedState.first, currentEstimatedState.second );

                    // Feed maneuver to the control system
                    controlSystem_->updateOrbitController( guidanceSystem_->getScheduledApsisManeuver( ) );
                }
                else
                {
                    // Check whether aerobraking is complete
                    isPropagationToBeStopped = isAerobrakingComplete( true );
                }

                // Run house keeping routines
                if ( navigationSystem_->currentOrbitCounter_ != 0 )
                {
                    runHouseKeepingRoutines( );
                }

                // Step up orbit counter
                navigationSystem_->currentOrbitCounter_++;
                guidanceSystem_->setCurrentOrbitCounter( navigationSystem_->currentOrbitCounter_ );

                // Renew random coefficients for perturbed atmosphere
                std::cerr << "Atmosphere randomization is OFF." << std::endl;
//                instrumentsModel_->randomizeAtmospherePerturbations( );

                // Invert completion flags
                maneuveringPhaseComplete_ = true;
                atmosphericPhaseComplete_ = false;

                // Inform user
                std::string orbitNumber = std::to_string( navigationSystem_->currentOrbitCounter_ - 1 );
                std::cout << std::endl << "-------------- ORBIT " << orbitNumber << " COMPLETED --------------" << std::endl;
            }
            // Check aerobraking phase and true anomaly to see if periapsis maneuvering phase
            else if ( guidanceSystem_->getIsAerobrakingPhaseActive( GuidanceSystem::termination_phase ) &&
                      ( currentEstimatedTrueAnomaly < ( 0.95 * PI ) ) && !atmosphericPhaseComplete_ )
            {
                // Inform user
                std::cout << std::endl << "REACHED FINAL PERIAPSIS. Preparing to perform maneuver." << std::endl;

                // Run periapsis maneuver estimator
                guidanceSystem_->runPeriapsisManeuverEstimator( currentTime, currentEstimatedState.first, currentEstimatedState.second );

                // Feed maneuver to the control system
                controlSystem_->updateOrbitController( guidanceSystem_->getScheduledApsisManeuver( ) );

                // Propagation is stopped to perform periapsis maneuver
                isPropagationToBeStopped = true;

                // Invert completion flags
                maneuveringPhaseComplete_ = false;
                atmosphericPhaseComplete_ = true;
            }
            // Check DAIA and true anomaly to see if post-atmosphere processes phase
            else if ( ( navigationSystem_->getIsSpacecraftAboveDynamicAtmosphericInterfaceAltitude( ) &&
                        ( currentEstimatedTrueAnomaly < ( 0.95 * PI ) ) ) && !atmosphericPhaseComplete_ )
            {
                // Inform user
                std::cout << std::endl << "EXITED ATMOSPHERE. Running post-atmosphere processes." << std::endl;

                // Retireve history of inertial measurement unit measurements
                std::map< double, Eigen::Vector6d > currentOrbitHistoryOfInertialMeasurementUnitMeasurements =
                        instrumentsModel_->getCurrentOrbitHistoryOfInertialMeasurmentUnitMeasurements( );

                // Extract measured translational accelerations and transform to inertial frame
                std::map< double, Eigen::Vector3d > currentOrbitHistoryOfMeasuredTranslationalAccelerations;
                for ( measurementConstantIterator_ = currentOrbitHistoryOfInertialMeasurementUnitMeasurements.begin( );
                      measurementConstantIterator_ != currentOrbitHistoryOfInertialMeasurementUnitMeasurements.end( );
                      measurementConstantIterator_++ )
                {
                    currentOrbitHistoryOfMeasuredTranslationalAccelerations[ measurementConstantIterator_->first ] =
                            measurementConstantIterator_->second.segment( 0, 3 );
                }

                // Perform periapse time and atmosphere estimations
                if ( propagators::IMAN_ANALYSIS_INDEX != 0 )
                {
                    std::cerr << "Post-atmosphere processes are OFF." << std::endl;
                }
                else
                {
                    navigationSystem_->runPostAtmosphereProcesses( currentOrbitHistoryOfMeasuredTranslationalAccelerations );
                }

                // Invert completion flags
                maneuveringPhaseComplete_ = false;
                atmosphericPhaseComplete_ = true;
            }

            // Save current value of propagation termination index
            previousIsPropagationToBeStopped_ = isPropagationToBeStopped;
            navigationSystem_->setCurrentTime( currentTime );
        }
        else
        {
            // Inform user
            if ( propagators::IMAN_ANALYSIS_INDEX != 1 )
            {
                std::cerr << "Warning in onboard computer. The current time (" << currentTime - initialTime_ << ") has already been " <<
                             "processed. Navigation time: " << navigationSystem_->getCurrentTime( ) - initialTime_ << "." << std::endl
                          << "Time difference: " << currentTime - navigationSystem_->getCurrentTime( ) << "." << std::endl;
            }

            // Return previous value of propagation termination index
            isPropagationToBeStopped = previousIsPropagationToBeStopped_;
            previousIsPropagationToBeStopped_ = false;
        }

        // Give output
        return isPropagationToBeStopped;
    }

    //! Function to check whether the aerobraking maneuver has been completed.
    /*!
     *  Function to check whether the aerobraking maneuver has been completed.
     *  \param isCallInternal Boolean denoting whether the function is being called internally; default value is false; the only
     *      difference is that if the call is internal, no messages will be printed.
     *  \return Boolean denoting whether the aerobraking maneuver has come to an end.
     */
    bool isAerobrakingComplete( const bool isCallInternal = false )
    {
        // Check if aerobraking is complete
        bool aerobrakingComplete = guidanceSystem_->getIsAerobrakingPhaseActive( GuidanceSystem::aerobraking_complete );

        // Inform user
        if ( aerobrakingComplete && !isCallInternal )
        {
            std::string velocityChange = std::to_string( guidanceSystem_->getHistoryOfApsisManeuverMagnitudes( ).first );
            std::cout << std::endl << "~~~~~~~~~~~~~~ AEROBRAKING COMPLETE ~~~~~~~~~~~~~~" << std::endl << std::endl
                      << "Cumulative change in velocity: " << velocityChange << " m/s" << std::endl;
        }

        // Give output
        return aerobrakingComplete;
    }

private:

    //! Function to run house keeping routines when new orbit is initiated.
    /*!
     *  Function to run house keeping routines when new orbit is initiated. This function erases the history of all variables that
     *  are only stored for the orbit. In this way, the computer frees-up space for other data.
     */
    void runHouseKeepingRoutines( )
    {
        // Empty maps and vectors of data belonging to the current orbit
        controlSystem_->clearCurrentOrbitControlHistory( );
        instrumentsModel_->clearCurrentOrbitMeasurementHistories( );
        navigationSystem_->clearCurrentOrbitEstimationHistory( );
    }

    //! Function to revert to the previous time step.
    /*!
     *  Function to revert to the previous time step. This function is run if the current propagation needs to be stopped, since
     *  the current time will be run the next time the GNC system is called.
     *  \param currentTime Double denoting the current time, i.e., the instant that has to be discarded.
     */
    void revertToPreviousTimeStep( const double currentTime )
    {
        // Empty elements belonging to the current time
        controlSystem_->revertToPreviousTimeStep( currentTime );
        instrumentsModel_->revertToPreviousTimeStep( currentTime );
        navigationSystem_->revertToPreviousTimeStep( currentTime );
    }

    //! Pointer to the control system for the aerobraking maneuver.
    const boost::shared_ptr< ControlSystem > controlSystem_;

    //! Pointer to the guidance system for the aerobraking maneuver.
    const boost::shared_ptr< GuidanceSystem > guidanceSystem_;

    //! Pointer to the navigation system for the aerobraking maneuver.
    const boost::shared_ptr< NavigationSystem > navigationSystem_;

    //! Pointer to the navigation system for the aerobraking maneuver.
    const boost::shared_ptr< InstrumentsModel > instrumentsModel_;

    //! Boolean denoting whether the maneuvering phase for this orbit has been complete.
    bool maneuveringPhaseComplete_;

    //! Boolean denoting whether the atmospheric phase for this orbit has been complete.
    bool atmosphericPhaseComplete_;

    //! Boolean denoting whether the maneuver should be performed, by adding the velocity change to the current state.
    bool performManeuverOnNextCall_;

    //! Pair denoting whether the Deep Space Network tracking is to be performed and the last day it was performed.
    std::pair< bool, unsigned int > deepSpaceNetworkTrackingInformation_;

    //! Double denoting initial time of the simulation.
    double initialTime_;

    //! Boolean denoting the previous flag for propagation stop.
    bool previousIsPropagationToBeStopped_;

    //! Predefined iterator to save (de)allocation time.
    std::map< double, Eigen::Vector6d >::const_iterator measurementConstantIterator_;

};

} // namespace thesis

} // namespace tudat

#endif // TUDAT_ONBOARD_COMPUTER_MODEL_H
