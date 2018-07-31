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

#include <iostream>
#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"

//! Typedefs and using statements to simplify code.
namespace Eigen { typedef Eigen::Matrix< double, 12, 1 > Vector12d; }

namespace tudat
{

namespace guidance_navigation_control
{

//! Class for control system of an aerobraking maneuver.
class ControlSystem
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param proportionalGain
     *  \param integralGain
     *  \param derivativeGain
     */
    ControlSystem( const Eigen::Vector3d& proportionalGain, const Eigen::Vector3d& integralGain,
                   const Eigen::Vector3d& derivativeGain ) :
        proportionalGain_( proportionalGain ), integralGain_( integralGain ), derivativeGain_( derivativeGain )
    { }

    //! Destructor.
    ~ControlSystem( ) { }

    //! Function to update the orbit controller with the scheduled apoapsis maneuver, computed by the guidance system.
    void updateOrbitController( const Eigen::Vector3d& scheduledApsoapsisManeuver )
    {
        // Set apoapsis maneuver magnitude and direction
        scheduledApsoapsisManeuver_ = scheduledApsoapsisManeuver;
    }

    //! Function to retireve current control vector for attitude.
    Eigen::Vector3d getCurrentAttitudeControlVector( ) { return Eigen::Vector3d::Zero( ); }

    //! Function to retirieve the apoapsis maneuver.
    Eigen::Vector3d getScheduledApoapsisManeuver( ) { return scheduledApsoapsisManeuver_; }

    //! Clear history of control vectors for the current orbit.
    void clearCurrentOrbitControlHistory( ) { }

private:

    //! Double denoting the proportional gain for the PID attitude controller.
    const Eigen::Vector3d proportionalGain_;

    //! Double denoting the integral gain for the PID attitude controller.
    const Eigen::Vector3d integralGain_;

    //! Double denoting the derivative gain for the PID attitude controller.
    const Eigen::Vector3d derivativeGain_;

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsoapsisManeuver_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_CONTROL_SYSTEM_H