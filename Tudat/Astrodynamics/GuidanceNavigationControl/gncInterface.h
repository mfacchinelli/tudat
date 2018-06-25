/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef MICHELE_GNC_INTERFACE
#define MICHELE_GNC_INTERFACE

#include <map>
#include <eigen/Core>
#include <boost/function.hpp>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/onboardComputer.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Class for the onboard computer of the spacecraft.
class GNCInterface
{
public:

    //! Constructor.
    GNCInterface( ) { }

    //! Destructor.
    virtual ~GNCInterface( ) { }

protected:

};

} // namespace thesis

} // namespace tudat

#endif // MICHELE_GNC_INTERFACE
