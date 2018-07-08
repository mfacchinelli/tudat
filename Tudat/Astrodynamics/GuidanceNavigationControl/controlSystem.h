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
namespace Eigen
{
typedef Eigen::Matrix< double, 16, 1 > Vector16d;
}

namespace tudat
{

namespace guidance_navigation_control
{

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
Eigen::Vector4d computeErrorInEstimatedQuaternion( const Eigen::Vector4d& currentQuaternionToBaseFrame,
                                                   const Eigen::Vector4d& currentCommandedQuaternionToBaseFrame );

//! Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
/*!
 * Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 * \param currentQuaternionsToBaseFrame Quaternions (in vector representation) that defined the rotation from body-fixed to inertial
 * frame.
 * \param angularVelocityVectorInBodyFixedFrame Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 */
Eigen::Vector4d calculateQuaternionDerivative( const Eigen::Vector4d& currentQuaternionsToBaseFrame,
                                               const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

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
    {
        // Initialize variables to default values
        currentControlVector_ = Eigen::Vector3d::Zero( );
        previousCommandedQuaternionState_ = Eigen::Vector4d::Zero( );
        isSignToBeFlipped_ = false;
    }

    //! Destructor.
    ~ControlSystem( ) { }

    //! Function to update the attitude controller.
    /*!
     *  Function to update the attitude controller. The controller is based on the PID (proportional, integral, derivative)
     *  principle. Thus, first the error w.r.t. the commanded quaternion is computed, and then this value is directly used to
     *  determine the proportional term of the control. The integral term is found by integrating the error over time, and the
     *  derivative term is found by finding the error betwee the derivative of the current quaternion w.r.t. its commanded value.
     *  The commanded value of the derivative is zero. Note that the integral term is found by integrating numerically, with the
     *  extended Simpson formula. However, whenever the number of available points is less than or equal to 6, this function
     *  reverts to the simple trapezoidal integration.
     *  \param currentEstimatedStateVector Current estimated state vector.
     *  \param currentMeasuredRotationalVelocityVector Current measured rotational velocity vector.
     *  \param navigationRefreshStepSize Refresh step size of the navigation system.
     *  \param currentMeanMotion Current mean motion as provided by the navigation system.
     */
    void updateAttitudeController( const Eigen::Vector16d& currentEstimatedStateVector,
                                   const Eigen::Vector3d& currentMeasuredRotationalVelocityVector,
                                   const double navigationRefreshStepSize,
                                   const double currentMeanMotion )
    {
        // Determine commanded quaternion
        Eigen::Vector4d currentCommandedQuaternionState = ( isSignToBeFlipped_ ? -1.0 : 1.0 ) *
                computeCurrentCommandedQuaternionState( currentEstimatedStateVector );

        // Check if current quaternion signs match the previous time otherwise flip the sign
        bool doesSignHistoryMatch = currentCommandedQuaternionState.isApprox( -previousCommandedQuaternionState_, 1e-1 );
        if ( doesSignHistoryMatch )
        {
            currentCommandedQuaternionState *= -1.0;
        }
        isSignToBeFlipped_ = isSignToBeFlipped_ || doesSignHistoryMatch;
        previousCommandedQuaternionState_ = currentCommandedQuaternionState;

        // Compute difference between current and commanded state
        Eigen::Vector4d currentErrorInEstimatedQuaternionState =
                computeErrorInEstimatedQuaternion( currentEstimatedStateVector.segment( 6, 4 ),
                                                   currentCommandedQuaternionState );
        historyOfQuaternionStateErrors_.push_back( currentErrorInEstimatedQuaternionState.segment( 1, 3 ) );

        // Compute difference between current and commanded derivative
        Eigen::Vector4d currentCommandedQuaternionDerivative = computeCurrentCommandedQuaternionDerivative( currentEstimatedStateVector,
                                                                                                            currentCommandedQuaternionState,
                                                                                                            currentMeanMotion );
        Eigen::Vector4d currentErrorInEstimatedQuaternionDerivative =
                computeErrorInEstimatedQuaternion( currentEstimatedStateVector.segment( 6, 4 ),
                                                   currentCommandedQuaternionDerivative ) +
                computeErrorInEstimatedQuaternion( calculateQuaternionDerivative( currentEstimatedStateVector.segment( 6, 4 ),
                                                                                  currentMeasuredRotationalVelocityVector ),
                                                   currentCommandedQuaternionState );

//        std::cout << "Current estimated quaternion: " << currentEstimatedStateVector.segment( 6, 4 ).transpose( ) << std::endl
//                  << "Current estimated derivative: "  <<
//                     calculateQuaternionDerivative( currentEstimatedStateVector.segment( 6, 4 ),
//                                                    currentMeasuredRotationalVelocityVector ).transpose( ) << std::endl
//                  << "Commanded state: " << currentCommandedQuaternionState.transpose( ) << std::endl
//                  << "Commnaded derivative: " << currentCommandedQuaternionDerivative.transpose( ) << std::endl
//                  << "Proportional: " <<
//                     proportionalGain_.cwiseProduct( currentErrorInEstimatedQuaternionState.segment( 1, 3 ) ).transpose( ) << std::endl
//                  << "Integral: " <<
//                     integralGain_.cwiseProduct( numerical_quadrature::performExtendedSimpsonsQuadrature(
//                                                     navigationRefreshStepSize, historyOfQuaternionStateErrors_ ) ).transpose( ) << std::endl
//                  << "Derivative: " <<
//                     derivativeGain_.cwiseProduct( currentErrorInEstimatedQuaternionDerivative.segment( 1, 3 ) ).transpose( ) << std::endl;

        // Compute control vector based on control gains and error
        currentControlVector_ = - ( proportionalGain_.cwiseProduct( currentErrorInEstimatedQuaternionState.segment( 1, 3 ) ) +
                                    integralGain_.cwiseProduct( numerical_quadrature::performExtendedSimpsonsQuadrature(
                                                                    navigationRefreshStepSize, historyOfQuaternionStateErrors_ ) ) +
                                    derivativeGain_.cwiseProduct( currentErrorInEstimatedQuaternionDerivative.segment( 1, 3 ) ) );
        currentOrbitHistoryOfControlVectors_.push_back( currentControlVector_ );
        // only the imaginary part of the quaternion is used, since only three terms are needed to fully control the spacecraft
//        std::cout << "Current control vector: " << currentControlVector_.transpose( ) << std::endl << std::endl;
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

    //! Function to retrieve history of control vectors for the current orbit.
    std::vector< Eigen::Vector3d > getCurrentOrbitHistoryOfControlVectors( ) { return currentOrbitHistoryOfControlVectors_; }

    //! Clear history of control vectors for the current orbit.
    void clearCurrentOrbitControlHistory( )
    {
        currentOrbitHistoryOfControlVectors_.clear( );
    }

private:

    //! Function to compute the current commanded quaternion to base frame.
    /*!
     *  Function to compute the current commanded quaternion to base frame. The body-fixed frame is assumed to correspond
     *  to the trajectory frame. Thus, the direction cosine matrix (DCM) describing the rotation from trajectory to inertial
     *  frame can be found and taken as full rotation from body-fixed to inertial. The DCM is found by using the velocity and
     *  radial distance vector. The velocity (unit) vector corresponds directly to the x-axis of the trajectory frame, whereas
     *  the z-axis is computed by subtracting from the radial distance (unit) vector its projection on the x-axis. Then, the
     *  y-axis is determined via the right-hand rule (i.e., with the cross product). The commanded state thus corresponds to a
     *  state with zero angle of attack, angle of side-slip and bank angle. Note that since the estimated state is used, the actual
     *  transformation can differ.
     *  \param currentEstimatedStateVector Current estimated state as provided by the navigation system.
     *  \return Quaternion representing the estimated rotation from trajectory to inertial frame. Thus the commanded
     *      quaternion corresponds to a state with zero angle of attack, angle of side-slip and bank angle.
     */
    Eigen::Vector4d computeCurrentCommandedQuaternionState( const Eigen::Vector16d& currentEstimatedStateVector )
    {
        // Declare direction cosine matrix
        Eigen::Matrix3d transformationFromInertialToTrajectoryFrame;

        // Find the trajectory x-axis unit vector
        Eigen::Vector3d xUnitVector = currentEstimatedStateVector.segment( 3, 3 ).normalized( );
        transformationFromInertialToTrajectoryFrame.col( 0 ) = xUnitVector;

        // Find trajectory z-axis unit vector
        Eigen::Vector3d zUnitVector = currentEstimatedStateVector.segment( 0, 3 ).normalized( );
        zUnitVector -= zUnitVector.dot( xUnitVector ) * xUnitVector;
        transformationFromInertialToTrajectoryFrame.col( 2 ) = zUnitVector;

        // Find body-fixed y-axis unit vector
        transformationFromInertialToTrajectoryFrame.col( 1 ) = zUnitVector.cross( xUnitVector );

        // Transform DCM to quaternion and give output
        return linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( transformationFromInertialToTrajectoryFrame ) );
        // note that due to the different rotation convention in Eigen, transforming the DCM with the constructor above
        // automatically takes the inverse of the rotation
    }

    //! Function to compute the current commanded quaternion derivative to base frame.
    /*!
     *  Function to compute the current commanded quaternion derivative to base frame. The derivative is computed by
     *  assuming that the expected rotational velocity of the spacecraft is equal to the mean motion of the spacecraft.
     *  With this assumption, the spacecraft attitude does not exactly match the commanded attitude (for elliptical orbits),
     *  but it gives a close approximation, especially for the regions near peri- and apoapsis.
     *  \param currentEstimatedStateVector Current estimated state as provided by the navigation system.
     *  \param transformationFromTrajectoryToInertialFrame
     *  \param currentMeanMotion Current mean motion as provided by the navigation system.
     *  \return Quaternion derivative representing the rotational rate equal to the mean motion of the spacecraft, applied
     *      around the y-axis (body-fixed).
     */
    Eigen::Vector4d computeCurrentCommandedQuaternionDerivative( const Eigen::Vector16d& currentEstimatedStateVector,
                                                                 const Eigen::Vector4d& transformationFromInertialToTrajectoryFrame,
                                                                 const double currentMeanMotion )
    {
        // Compute the mean motion in trajectory frame
        Eigen::Vector3d expectedRotationalVelocityVector = Eigen::Vector3d::Zero( );
        expectedRotationalVelocityVector[ 1 ] = currentMeanMotion;

        // Transform the rotational velocity vector
        expectedRotationalVelocityVector = (
                    linear_algebra::convertVectorToQuaternionFormat( currentEstimatedStateVector.segment( 6, 4 ) ).toRotationMatrix( ) *
                    linear_algebra::convertVectorToQuaternionFormat( transformationFromInertialToTrajectoryFrame ).toRotationMatrix( ) ) *
                expectedRotationalVelocityVector;

        // Compute the expected derivative
        return calculateQuaternionDerivative( currentEstimatedStateVector.segment( 6, 4 ), expectedRotationalVelocityVector );
    }

    //! Double denoting the proportional gain for the PID attitude controller.
    const Eigen::Vector3d proportionalGain_;

    //! Double denoting the integral gain for the PID attitude controller.
    const Eigen::Vector3d integralGain_;

    //! Double denoting the derivative gain for the PID attitude controller.
    const Eigen::Vector3d derivativeGain_;

    //! Vector denoting the previous signs of quaternions.
    Eigen::Vector4d previousCommandedQuaternionState_;

    //! Boolean denoting whether the sign of the quaternion needs to be flipped.
    bool isSignToBeFlipped_;

    //! Vector denoting the current quaternion attitude correction.
    Eigen::Vector3d currentControlVector_;

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsoapsisManeuver_;

    //! History of errors in the estimated quaternion state.
    std::vector< Eigen::Vector3d > historyOfQuaternionStateErrors_;

    //! History of control torque vectors for current orbit.
    std::vector< Eigen::Vector3d > currentOrbitHistoryOfControlVectors_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_CONTROL_SYSTEM_H
