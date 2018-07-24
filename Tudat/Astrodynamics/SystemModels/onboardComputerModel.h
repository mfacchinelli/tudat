/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ONBOARD_COMPUTER_MODEL_H
#define TUDAT_ONBOARD_COMPUTER_MODEL_H

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/GuidanceNavigationControl/controlSystem.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/guidanceSystem.h"

#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"

//! Typedefs and using statements to simplify code.
namespace Eigen { typedef Eigen::Matrix< double, 16, 1 > Vector16d; }
using namespace tudat::guidance_navigation_control;

namespace tudat
{

namespace system_models
{

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector16d onboardSystemModel( const double currentTime, const Eigen::Vector16d& currentEstimatedStateVector,
                                     const Eigen::Vector3d& currentControlVector,
                                     const boost::function< Eigen::Vector3d( const Eigen::Vector6d& ) >&
                                     currentEstimatedTranslationalAccelerationFunction,
                                     const Eigen::Vector3d& currentMeasuredRotationalVelocityVector );

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector7d onboardMeasurementModel( const double currentTime, const Eigen::Vector16d& currentEstimatedStateVector,
                                         const boost::function< Eigen::Vector3d( const Eigen::Vector6d& ) >&
                                         currentEstimatedNonGravitationalTranslationalAccelerationFunction );

//! Class for the onboard computer of the spacecraft.
class OnboardComputerModel
{
public:

    //! Constructor.
    OnboardComputerModel( const boost::shared_ptr< ControlSystem > controlSystem,
                          const boost::shared_ptr< GuidanceSystem > guidanceSystem,
                          const boost::shared_ptr< NavigationSystem > navigationSystem,
                          const boost::shared_ptr< NavigationInstrumentsModel > instrumentsModel ) :
        controlSystem_( controlSystem ), guidanceSystem_( guidanceSystem ), navigationSystem_( navigationSystem ),
        instrumentsModel_( instrumentsModel )
    {
        dummyCallCounter_ = 0; // <<<<<<<<<----------

        // Define internal constants
        maneuveringPhaseComplete_ = true; // simulation starts at apoapsis with no need to perform a maneuver
        atmosphericPhaseComplete_ = false;
        atmosphericInterfaceRadius_ = navigationSystem_->getAtmosphericInterfaceRadius( );

        // Create functions to retrieve accelerations from navigation system
        boost::function< Eigen::Vector3d( const Eigen::Vector6d& ) > translationalAccelerationFunction =
                boost::bind( &NavigationSystem::getCurrentEstimatedTranslationalAcceleration, navigationSystem_, _1 );
        boost::function< Eigen::Vector3d( const Eigen::Vector6d& ) > nonGravitationalTranslationalAccelerationFunction =
                boost::bind( &NavigationSystem::getCurrentEstimatedNonGravitationalTranslationalAcceleration,
                             navigationSystem_, _1 );

        // Create navigation system objects
        navigationSystem_->createNavigationSystemObjects(
                    boost::bind( &onboardSystemModel, _1, _2,
                                 boost::bind( &ControlSystem::getCurrentAttitudeControlVector, controlSystem_ ),
                                 translationalAccelerationFunction,
                                 boost::bind( &NavigationInstrumentsModel::getCurrentGyroscopeMeasurement, instrumentsModel_ ) ),
                    boost::bind( &onboardMeasurementModel, _1, _2, nonGravitationalTranslationalAccelerationFunction ) );

        // Create guidance system objects
        guidanceSystem_->createGuidanceSystemObjects( boost::bind( &NavigationSystem::propagateStateWithCustomTerminationSettings,
                                                                        navigationSystem_, _1, _2, -1.0 ) );
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

        // Update instrument models and extract measurements
        instrumentsModel_->updateInstruments( currentTime );
        Eigen::Vector7d currentExternalMeasurementVector;
        currentExternalMeasurementVector.segment( 0, 3 ) = instrumentsModel_->getCurrentAccelerometerMeasurement( );
        currentExternalMeasurementVector.segment( 3, 4 ) = instrumentsModel_->getCurrentStarTrackerMeasurement( );

        // Update filter from previous time to next time
        navigationSystem_->runStateEstimator( currentTime, currentExternalMeasurementVector,
                                              instrumentsModel_->getCurrentGyroscopeMeasurement( ) );

        // Update attitude controller
        controlSystem_->updateAttitudeController(
                    navigationSystem_->getCurrentEstimatedState( ),
                    removeErrorsFromInertialMeasurementUnitMeasurement( instrumentsModel_->getCurrentGyroscopeMeasurement( ),
                                                                        navigationSystem_->getCurrentEstimatedState( ).segment( 10, 6 ) ),
                    navigationSystem_->getNavigationRefreshStepSize( ),
                    navigationSystem_->getCurrentEstimatedMeanMotion( ) );

        // Check if stopping condition is met or if the post-atmospheric phase processes need to be carried out
        std::pair< Eigen::VectorXd, Eigen::VectorXd > currentEstimatedState = navigationSystem_->getCurrentEstimatedTranslationalState( );
        double currentEstimatedTrueAnomaly = currentEstimatedState.second[ 5 ];
        if ( ( currentEstimatedTrueAnomaly >= PI ) && !maneuveringPhaseComplete_ ) // check true anomaly
        {
            // Inform user
            std::cout << "Reached apoapsis. Preparing to perform maneuver." << std::endl;

            // Determine in which phase of aerobraking the spacecraft is and perform corridor estimation
            guidanceSystem_->determineAerobrakingPhase( currentEstimatedState.second,
                                                        navigationSystem_->getAtmosphereInitiationIndicators( ) );
            guidanceSystem_->runCorridorEstimator( currentTime, currentEstimatedState.second,
                                                   navigationSystem_->getRadius( ),
                                                   navigationSystem_->getStandardGravitationalParameter( ) );

            // Check whether propagation needs to be stopped
            // Propagation is stopped only if an apoapsis maneuver needs to be performed
            isPropagationToBeStopped = guidanceSystem_->getIsApoapsisManeuverToBePerformed( );

            // If propagation is to be stopped, to add Delta V to state, run remaining pre-maneuver processes
            if ( isPropagationToBeStopped )
            {
                // Run maneuver estimator
                guidanceSystem_->runManeuverEstimator( currentEstimatedState.first, currentEstimatedState.second,
                                                       navigationSystem_->getCurrentEstimatedMeanMotion( ),
                                                       navigationSystem_->getRadius( ) );

                // Retireve required apoapsis maneuver
                Eigen::Vector3d scheduledApoapsisManeuver = guidanceSystem_->getScheduledApoapsisManeuver( );

                // Feed maneuver to the navigation system
                Eigen::Vector6d updatedCurrentEstimatedCartesianState = navigationSystem_->getCurrentEstimatedTranslationalState( ).first;
                updatedCurrentEstimatedCartesianState.segment( 3, 3 ) += scheduledApoapsisManeuver;
                navigationSystem_->setCurrentEstimatedCartesianState( updatedCurrentEstimatedCartesianState );

                // Feed maneuver to the control system
                controlSystem_->updateOrbitController( scheduledApoapsisManeuver );
            }

            // Run house keeping routines
            runHouseKeepingRoutines( );

            // Step up orbit counter
            navigationSystem_->currentOrbitCounter_++;

            // Invert completion flags
            maneuveringPhaseComplete_ = true;
            atmosphericPhaseComplete_ = false;
        }
        else if ( ( ( ( currentEstimatedState.first.segment( 0, 3 ).norm( ) - atmosphericInterfaceRadius_ ) > 0.0 ) &&
                    ( ( currentEstimatedTrueAnomaly >= 0.0 ) && ( currentEstimatedTrueAnomaly < ( 0.5 * PI ) ) ) ) &&
                  !atmosphericPhaseComplete_ ) // check altitude
        {
            // Inform user
            std::cout << "Exited atmosphere. Running post-atmosphere processes." << std::endl;

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
            navigationSystem_->runPostAtmosphereProcesses( currentOrbitHistoryOfMeasuredTranslationalAccelerations );

            // Invert completion flags
            maneuveringPhaseComplete_ = false;
            atmosphericPhaseComplete_ = true;
        }

        // Give output
        return isPropagationToBeStopped;
    }

    //! Function to check whether the aerobraking maneuver has been completed.
    /*!
     *  Function to check whether the aerobraking maneuver has been completed.
     *  \return Boolean denoting whether the aerobraking maneuver has been completed.
     */
    bool isAerobrakingComplete( )
    {
        dummyCallCounter_++;
        return ( dummyCallCounter_ > 0 );
    }

private:

    unsigned int dummyCallCounter_;

    //! Function to run house keeping routines when new orbit is initiated.
    void runHouseKeepingRoutines( )
    {
        // Empty maps and vectors of data belonging to the current orbit
        controlSystem_->clearCurrentOrbitControlHistory( );
        instrumentsModel_->clearCurrentOrbitMeasurementHistories( );
        navigationSystem_->clearCurrentOrbitEstimationHistory( );
    }

    //! Boolean denoting whether the maneuvering phase for this orbit has been complete.
    bool maneuveringPhaseComplete_;

    //! Boolean denoting whether the atmospheric phase for this orbit has been complete.
    bool atmosphericPhaseComplete_;

    //! Double denoting the atmospheric interface altitude.
    double atmosphericInterfaceRadius_;

    //! Pointer to the control system for the aerobraking maneuver.
    const boost::shared_ptr< ControlSystem > controlSystem_;

    //! Pointer to the guidance system for the aerobraking maneuver.
    const boost::shared_ptr< GuidanceSystem > guidanceSystem_;

    //! Pointer to the navigation system for the aerobraking maneuver.
    const boost::shared_ptr< NavigationSystem > navigationSystem_;

    //! Pointer to the navigation system for the aerobraking maneuver.
    const boost::shared_ptr< NavigationInstrumentsModel > instrumentsModel_;

    //! Predefined iterator to save (de)allocation time.
    std::map< double, Eigen::Vector6d >::const_iterator measurementConstantIterator_;

};

} // namespace thesis

} // namespace tudat

#endif // TUDAT_ONBOARD_COMPUTER_MODEL_H
