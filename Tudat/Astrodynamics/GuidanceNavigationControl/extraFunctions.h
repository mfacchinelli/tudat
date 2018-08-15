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
#include <tuple>
#include <iostream>

#include <boost/function.hpp>

#include "Tudat/Basics/basicTypedefs.h"

//! Typedefs and using statements to simplify code.
namespace Eigen { typedef Eigen::Matrix< double, 12, 1 > Vector12d; typedef Eigen::Matrix< double, 5, 5 > Matrix5d;
                  typedef Eigen::Matrix< double, 12, 12 > Matrix12d; }

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to compute the Jacobian matrix for the system function.
/*!
 *  Function to compute the Jacobian matrix for the system function. The accelerations considered are the gravitational spherical
 *  harmonics up to term (2,0) and the aerodynamic acceleration.
 *  \param currentTime Double denoting current time.
 *  \param currentState Vector denoting current state as given by filter.
 *  \param densityFunction Function returning the density corresponding to the input Cartesian state.
 *  \param planetGravitationalParameter Double denoting the gravitational parameter of the planet being orbited.
 *  \param aerodynamicParameter Double denoting the aerodynamic parameter of the spacecraft.
 *  \return Jacobian matrix for the system function at the current state.
 */
Eigen::Matrix12d computeSystemJacobianMatrix( const double currentTime, const Eigen::Vector12d& currentState,
                                              const boost::function< double( const Eigen::Vector6d& ) >& densityFunction,
                                              const double planetGravitationalParameter,
                                              const double aerodynamicParameter );

//! Function to compute the Jacobian matrix for the measurement function.
/*!
 *  Function to compute the Jacobian matrix for the measurement function.
 *  \param currentTime Double denoting current time.
 *  \param currentState Vector denoting current state as given by filter.
 *  \return Jacobian matrix for the measurement function at the current state.
 */
Eigen::Matrix3d computeMeasurementJacobianMatrix( const double currentTime, const Eigen::Vector12d& currentState );

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

//! Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
/*!
 *  Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
 *  \param currentAltitudeGuess
 *  \param initialEstimatedKeplerianState
 *  \param planetaryRadius
 *  \param planetaryGravitationalParameter
 *  \param maximumHeatRate
 *  \param maximumHeatLoad
 *  \param statePropagationFunction
 *  \return
 */
double lowerAltitudeBisectionFunction( const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
                                       const double planetaryRadius, const double planetaryGravitationalParameter,
                                       const double maximumHeatRate, const double maximumHeatLoad,
                                       const boost::function< std::pair< std::map< double, Eigen::VectorXd >,
                                       std::map< double, Eigen::VectorXd > >( const Eigen::Vector6d& ) >& statePropagationFunction );

//! Function to be used as input to the root-finder to determine the upper altitude bound for the periapsis corridor.
/*!
 *  Function to be used as input to the root-finder to determine the upper altitude bound for the periapsis corridor.
 *  \param currentAltitudeGuess
 *  \param initialEstimatedKeplerianState
 *  \param planetaryRadius
 *  \param planetaryGravitationalParameter
 *  \param minimumDynamicPressure
 *  \param statePropagationFunction
 *  \return
 */
double upperAltitudeBisectionFunction( const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
                                       const double planetaryRadius, const double planetaryGravitationalParameter,
                                       const double minimumDynamicPressure,
                                       const boost::function< std::pair< std::map< double, Eigen::VectorXd >,
                                       std::map< double, Eigen::VectorXd > >( const Eigen::Vector6d& ) >& statePropagationFunction );

//! Function to be used as input to the root-finder to determine the magnitude of the apoapsis maneuver.
/*!
 *  Function to be used as input to the root-finder to determine the magnitude of the apoapsis maneuver.
 *  \param currentMagnitudeGuess
 *  \param initialEstimatedCartesianState
 *  \param targetPeriapsisRadius
 *  \param transformationFromLocalToInertialFrame
 *  \param statePropagationFunction
 *  \return
 */
double maneuverBisectionFunction( const double currentMagnitudeGuess, const Eigen::Vector6d& initialEstimatedCartesianState,
                                  const double targetPeriapsisRadius, const Eigen::Matrix3d& transformationFromLocalToInertialFrame,
                                  const boost::function< std::pair< std::map< double, Eigen::VectorXd >,
                                  std::map< double, Eigen::VectorXd > >( const Eigen::Vector6d& ) >& statePropagationFunction );

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
