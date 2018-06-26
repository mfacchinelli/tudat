/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef MICHELE_GNC_COMPUTER
#define MICHELE_GNC_COMPUTER

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"

#include "Tudat/Astrodynamics/GuidanceNavigationControl/control.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/instrument.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigation.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/guidance.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Class for the onboard computer of the spacecraft.
class OnboardComputer
{
public:

    //! Constructor.
    OnboardComputer( const boost::shared_ptr< ControlSystem > controlSystem,
                     const boost::shared_ptr< GuidanceSystem > guidanceSystem,
                     const boost::shared_ptr< NavigationSystem > navigationSystem,
                     const boost::shared_ptr< InstrumentSystem > instrumentSystem ) :
        controlSystem_( controlSystem ), guidanceSystem_( guidanceSystem ), navigationSystem_( navigationSystem )
    {
        // Define internal constants
        maneuveringPhaseComplete_ = true;
        atmosphericPhaseComplete_ = false;
        atmosphericInterfaceAltitude_ = navigationSystem_->getRadius( ) + 500.0e3;
    }

    //! Destructor.
    ~OnboardComputer( ) { }

    //! Function to update the navigation system.
    /*!
     *  Function to check whether the propagation is to be be stopped, based on the estimated state. The propagation
     *  is stopped if either one of the two conditions below is met:
     *      - spacecraft is estimated to be at apoapsis;
     *      - spacecraft is estimated to have just left the atmosphere.
     *  To make sure that the propagation is checked every time step, this function needs to be fed to the propagator
     *  via the custom propagation termination settings.
     *  \param currentTime Current time in propagation.
     *  \return Boolean denoting whether the propagation has to be stopped.
     */
    bool checkStoppingCondition( const double currentTime )
    {
        // Update current time
        previousTime_ = navigationSystem_->getCurrentTime( );
        navigationSystem_->setCurrentTime( currentTime );

        // Update filter from previous time to next time
        navigationSystem_->stateEstimator( previousTime_ );

        // Check if either stopping condition is met
        std::pair< Eigen::VectorXd, Eigen::VectorXd > currentEstimatedState = navigationSystem_->getCurrentEstimatedState( );
        double currentEstimatedTrueAnomaly = orbital_element_conversions::convertEllipticalEccentricAnomalyToTrueAnomaly(
                    orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly( currentEstimatedState.second[ 2 ],
                    currentEstimatedState.second[ 5 ] ), currentEstimatedState.second[ 2 ] );
        if ( currentEstimatedTrueAnomaly >= mathematical_constants::PI && !maneuveringPhaseComplete_ ) // check mean anomaly
        {

            // Invert completion flags
            maneuveringPhaseComplete_ = true;
            atmosphericPhaseComplete_ = false;
        }
        else if ( ( ( currentEstimatedState.first.segment( 0, 3 ).norm( ) - atmosphericInterfaceAltitude_ ) > 0.0 &&
                    currentEstimatedTrueAnomaly >= 0.0 ) && !atmosphericPhaseComplete_ ) // check altitude
        {

            // Invert completion flags
            maneuveringPhaseComplete_ = false;
            atmosphericPhaseComplete_ = true;
        }
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
    boost::shared_ptr< ControlSystem > controlSystem_;

    //! Pointer to the guidance system for the aerobraking maneuver.
    boost::shared_ptr< GuidanceSystem > guidanceSystem_;

    //! Pointer to the navigation system for the aerobraking maneuver.
    boost::shared_ptr< NavigationSystem > navigationSystem_;

};

} // namespace thesis

} // namespace tudat

#endif // MICHELE_GNC_COMPUTER
