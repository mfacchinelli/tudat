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

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"

//! Typedefs and using statements to simplify code.
typedef Eigen::Matrix< double, 16, 1 > Eigen::Vector16d;

//! Function to compute the error in the current quaternion state, based on the commanded quaternion.
/*!
 *  Function to compute the error in the current quaternion state, based on the commanded quaternion. The
 *  error quaternion is computed by taking the difference between the commanded and estimated quaternions.
 *  Note that the real part of the error quaternion is zero by default, since only three of the quaternion
 *  elements are needed to compute the control vector.
 *  \param currentQuaternionToBaseFrame Current estimated quaternion, expressing the rotation from the
 *      body-fixed frame to the inertial frame.
 *  \param currentCommandedQuaternionToBaseFrame Current commanded quaternion, expressing the wanted rotation
 *      from body-fixed frame to the inertial frame.
 *  \return Error quaternion, expressing the difference between the commanded and the estimated quaternions.
 */
Eigen::Vector4d computeErrorInEstimatedQuaternionState( const Eigen::Vector4d& currentQuaternionToBaseFrame,
                                                        const Eigen::Vector4d& currentCommandedQuaternionToBaseFrame );

namespace tudat
{

namespace guidance_navigation_control
{

//! Class for control system of an aerobraking maneuver.
class ControlSystem
{
public:

    //! Constructor.
    ControlSystem( const double proportionalGain, const double integralGain, const double derivativeGain,
                   const double  ) :
        proportionalGain_( proportionalGain ), integralGain_( integralGain ), derivativeGain_( derivativeGain ),
        navigationRefreshStepSize_( TUDAT_NAN )
    {
        // Initialize control vector to zero
        currentControlVector_ = Eigen::Vector3d::Zero( );
    }

    //! Destructor.
    ~ControlSystem( ) { }

    //! Attitude control system.
    void updateAttitudeController( const Eigen::Vector16d& currentEstimatedState,
                                   const Eigen::Vector3d& currentMeasuredRotationalVelocityVector,
                                   const double navigationRefreshStepSize )
    {
        // Compute difference between current and commanded state
        Eigen::Vector4d currentErrorInEstimatedQuaternionState =
                computeErrorInEstimatedQuaternionState( currentEstimatedState.segment( 0, 4 ),
                                                        computeCurrentCommandedQuaternionState( currentEstimatedState ) );
        historyOfQuaternionStateErrors_.push_back( currentErrorInEstimatedQuaternionState );

        // Compute control vector based on control gains and error
        currentControlVector_ = proportionalGain_ * currentErrorInEstimatedQuaternionState +
                integralGain_ * numerical_quadrature::performExtendedSimpsonsQuadrature( navigationRefreshStepSize,
                                                                                         historyOfQuaternionStateErrors_ ) +
                derivativeGain_ * propagators::calculateQuaternionsDerivative( currentEstimatedState.segment( 0, 4 ),
                                                                               currentMeasuredRotationalVelocityVector );
    }

    //! Function to update the orbit controller with the scheduled apoapsis maneuver, computed by the guidance system.
    void updateOrbitController( const Eigen::Vector3d& scheduledApsoapsisManeuver )
    {
        // Set apoapsis maneuver magnitude and direction
        scheduledApsoapsisManeuver_ = scheduledApsoapsisManeuver;
    }

    //! Function to retireve current control vector for attitude.
    Eigen::Vector3d getCurrentAttitudeControlVector( ) { return currentControlVector_; }

    //! Function to retirieve the apoapsis maneuver.
    Eigen::Vector3d getScheduledApoapsisManeuver( ) { return scheduledApsoapsisManeuver_; }

private:

    //! Function to compute the current commanded quaternion to base frame.
    Eigen::Vector4d computeCurrentCommandedQuaternionState( const Eigen::Vector16d& currentEstimatedState )
    {
        return Eigen::Vector4d::UnitX( );
    }

    //! Double denoting the proportional gain for the PID attitude controller.
    const double proportionalGain_;

    //! Double denoting the integral gain for the PID attitude controller.
    const double integralGain_;

    //! Double denoting the derivative gain for the PID attitude controller.
    const double derivativeGain_;

    //! Vector denoting the current quaternion attitude correction.
    Eigen::Vector3d currentControlVector_;

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsoapsisManeuver_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_CONTROL_SYSTEM_H
