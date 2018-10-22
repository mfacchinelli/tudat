/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_IMAN_RMS_METHOD_H
#define TUDAT_IMAN_RMS_METHOD_H

#include <complex>
#include <cmath>
#include <limits>

namespace tudat
{

namespace propagators
{

// IMAN analysis index:
//      0: nominal simulation
//      1: RMS simulation
//      2: tuning (loop) simulation
extern unsigned int IMAN_ANALYSIS_INDEX;

} // namespace propagators

} // namespace tudat

#endif // TUDAT_IMAN_RMS_METHOD_H
