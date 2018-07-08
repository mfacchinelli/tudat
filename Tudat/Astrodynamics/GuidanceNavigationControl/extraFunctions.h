/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GNC_EXTRA_FUNCTIONS_H
#define TUDAT_GNC_EXTRA_FUNCTIONS_H

#include <map>
#include <iostream>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
/*!
 *  Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
 *  The function returns the difference between the areas under the acceleration curve, computed before and after
 *  a time guess provided by the root-finder object. Since the area under the curve is a constant value, the
 *  slicing procedure will always result in a larger and a smaller value of area, unless the slicing occurs exactly
 *  at the centroid of the area curve (which is the tagert point). Since no interpolation is used, it is possible that
 *  the time guess is not cointained in the onboardTime vector, thus, the nearest index is used as reference to
 *  compute the area.
 *  \param currentTimeGuess Current guess in periapse time to be used as slicing parameter.
 *  \param constantTimeStep Value of constant time step to be used for numerical integration.
 *  \param onboardTime Onboard times below the atmospheric interface.
 *  \param estimatedAerodynamicAccelerationMagnitude Vector of estimated aerodynamic acceleration magnitudes below the
 *      atmospheric interface altitude.
 *  \return Double representing the difference between the areas under the acceleration curve computed before and
 *      after the current time estimate (i.e., the slicing parameter).
 */
double areaBisectionFunction( const double currentTimeGuess, const double constantTimeStep,
                              const Eigen::VectorXd& onboardTime,
                              const std::vector< double >& estimatedAerodynamicAccelerationMagnitude );

//! Function to be used as input to the non-linear least squares process to determine the accelerometer errors.
/*!
 *  Function to be used as input to the non-linear least squares process to determine the accelerometer errors.
 *  \param currentErrorEstimate
 *  \param vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface Vector of estimated acceleration
 *  \return Pair of expected acceleration and Jacobian of the acceleration function w.r.t. the accelerometer errors.
 */
std::pair< Eigen::VectorXd, Eigen::MatrixXd > accelerometerErrorEstimationFunction(
        const Eigen::Vector6d& currentErrorEstimate,
        const std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface );

//! Function to be used as input to the non-linear least squares process to determine the atmosphere parameters of the
//! three-term atmosphere model.
std::pair< Eigen::VectorXd, Eigen::MatrixXd > threeModelParametersEstimationFunction(
        const Eigen::Vector5d& currentParameterEstimate,
        const Eigen::VectorXd& vectorOfEstimatedAltitudesBelowAtmosphericInterface,
        const double referenceAltitude );

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_GNC_EXTRA_FUNCTIONS_H
