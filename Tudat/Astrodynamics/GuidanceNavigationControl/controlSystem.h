/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONTROL_SYSTEM_H
#define TUDAT_CONTROL_SYSTEM_H

#include <Eigen/Core>

namespace tudat
{

namespace guidance_navigation_control
{

//! Class for control system of an aerobraking maneuver.
class ControlSystem
{
public:

    //! Constructor.
    ControlSystem( const double proportionalGain, const double integralGain, const double derivativeGain ) :
        proportionalGain_( proportionalGain ), integralGain_( integralGain ), derivativeGain_( derivativeGain )
    { }

    //! Destructor.
    ~ControlSystem( ) { }

    //! Attitude control system.
    void updateAttitudeController( )
    {
        currentControlVector_ = Eigen::Vector3d::Zero( );
    }

    //! Function to update the orbit controller with the scheduled apoapsis maneuver, computed by the guidance system.
    void updateOrbitController( const Eigen::Vector3d& scheduledApsoapsisManeuver )
    {
        scheduledApsoapsisManeuver_ = scheduledApsoapsisManeuver;
    }

    //! Function to retireve current control vector for attitude.
    Eigen::Vector3d getCurrentAttitudeControlVector( ) { return currentControlVector_; }

    //! Function to retirieve the apoapsis maneuver.
    Eigen::Vector3d getScheduledApoapsisManeuver( ) { return scheduledApsoapsisManeuver_; }

private:

    //! Vector denoting the current quaternion attitude correction.
    Eigen::Vector3d currentControlVector_;

    //! Double denoting the proportional gain for the PID attitude controller.
    const double proportionalGain_;

    //! Double denoting the integral gain for the PID attitude controller.
    const double integralGain_;

    //! Double denoting the derivative gain for the PID attitude controller.
    const double derivativeGain_;

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsoapsisManeuver_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_CONTROL_SYSTEM_H
