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

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/GuidanceNavigationControl/controlSystem.h"
#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/guidanceSystem.h"

namespace tudat
{

namespace guidance_navigation_control
{

using namespace tudat::system_models;

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
        // Define internal constants
        maneuveringPhaseComplete_ = true; // simulation starts at apoapsis with no need to perform a maneuver
        atmosphericPhaseComplete_ = false;
        atmosphericInterfaceAltitude_ = navigationSystem_->getRadius( ) + 500.0e3;
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
    bool checkStoppingCondition( const double currentTime )
    {
        // Define output value
        isPropagationToBeStopped = false;

        // Update current time
        previousTime_ = navigationSystem_->getCurrentTime( );
        navigationSystem_->setCurrentTime( currentTime );

        // Update filter from previous time to next time
        navigationSystem_->runStateEstimator( previousTime_, instrumentsModel_ );

        // Check if stopping condition is met or if the post-atmospheric phase processes need to be carried out
        std::pair< Eigen::VectorXd, Eigen::VectorXd > currentEstimatedState = navigationSystem_->getCurrentEstimatedState( );
        double currentEstimatedTrueAnomaly = currentEstimatedState.second[ orbital_element_conversions::trueAnomalyIndex ];
        if ( currentEstimatedTrueAnomaly >= mathematical_constants::PI && !maneuveringPhaseComplete_ ) // check true anomaly
        {
            // Stop propagation to add Delta V to actual state
            isPropagationToBeStopped = true;

            // Retireve required apoapsis maneuver and feed it to the control system
            controlSystem_->updateOrbitController( guidanceSystem_->getScheduledApoapsisManeuver( ) );

            // Invert completion flags
            maneuveringPhaseComplete_ = true;
            atmosphericPhaseComplete_ = false;
        }
        else if ( ( ( currentEstimatedState.first.segment( 0, 3 ).norm( ) - atmosphericInterfaceAltitude_ ) > 0.0 &&
                    currentEstimatedTrueAnomaly >= 0.0 ) && !atmosphericPhaseComplete_ ) // check altitude
        {
            // Retireve history of inertial measurement unit measurements
            std::map< double, Eigen::Vector6d > currentOrbitHistoryOfInertialMeasurementUnitMeasurements =
                    instrumentsModel_->getCurrentOrbitHistoryOfInertialMeasurmentUnitMeasurements( );

            // Extract measured translational accelerations
            std::map< double, Eigen::Vector3d > currentOrbitHistoryOfEstimatedAccelerations;
            for ( measurementConstantIterator_ measurementIterator = currentOrbitHistoryOfInertialMeasurementUnitMeasurements.begin( );
                  measurementIterator != currentOrbitHistoryOfInertialMeasurementUnitMeasurements.end( );
                  measurementIterator++ )
            {
                currentOrbitHistoryOfEstimatedAccelerations[ measurementIterator->first ] = measurementIterator->second.segment( 0, 3 );
            }

            // Perform periapse time estimation
            navigationSystem_->runPeriapseTimeEstimator( currentOrbitHistoryOfEstimatedAccelerations );

            // Perform atmosphere estimation
            navigationSystem_->runAtmosphereEstimator( currentOrbitHistoryOfEstimatedAccelerations );

            // Perform corridor estimation
            guidanceSystem_->runCorridorEstimator( );

            // Perform maneuver estimation
            guidanceSystem_->runManeuverEstimator( );

            // Invert completion flags
            maneuveringPhaseComplete_ = false;
            atmosphericPhaseComplete_ = true;
        }

        // Give output
        return isPropagationToBeStopped;
    }

private:

    //! Double denoting the previous time in the estimation process.
    double previousTime_;

    //! Boolean denoting whether the maneuvering phase for this orbit has been complete.
    bool maneuveringPhaseComplete_;

    //! Boolean denoting whether the atmospheric phase for this orbit has been complete.
    bool atmosphericPhaseComplete_;

    //! Double denoting the atmospheric interface altitude.
    double atmosphericInterfaceAltitude_;

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
