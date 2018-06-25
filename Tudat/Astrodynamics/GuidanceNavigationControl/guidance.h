/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef MICHELE_GNC_GUIDANCE
#define MICHELE_GNC_GUIDANCE

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

    //! Corridor estimator (CE).
    void corridorEstimator( );

    //! Maneuver estimator (ME).
    void maneuverEstimator( );

private:

};

} // namespace thesis

} // namespace tudat

#endif // MICHELE_GNC_GUIDANCE
