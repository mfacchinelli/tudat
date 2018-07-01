/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GUIDANCE_SYSTEM_H
#define TUDAT_GUIDANCE_SYSTEM_H

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Class for guidance system of an aerobraking maneuver.
class GuidanceSystem
{
public:

    //! Constructor.
    GuidanceSystem( ) { }

    //! Destructor.
    ~GuidanceSystem( ) { }

    //! Function to run corridor estimator (CE).
    void runCorridorEstimator( );

    //! Function to run maneuver estimator (ME).
    void runManeuverEstimator( )
    {
        scheduledApsoapsisManeuver_ = Eigen::Vector3d::Zero( );
    }

    //! Function to retirieve the apoapsis maneuver.
    Eigen::Vector3d getScheduledApoapsisManeuver( ) { return scheduledApsoapsisManeuver_; }

private:

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsoapsisManeuver_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_GUIDANCE_SYSTEM_H
