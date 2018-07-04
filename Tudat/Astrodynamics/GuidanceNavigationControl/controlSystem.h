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
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

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
Eigen::Vector4d computeErrorInEstimatedQuaternionState( const Eigen::Vector4d& currentQuaternionToBaseFrame,
                                                        const Eigen::Vector4d& currentCommandedQuaternionToBaseFrame );

//! Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
/*!
 * Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 * \param currentQuaternionsToBaseFrame Quaternions (in vector representation) that defined the rotation from body-fixed to inertial
 * frame.
 * \param angularVelocityVectorInBodyFixedFrame Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 */
Eigen::Vector4d calculateQuaternionsDerivative( const Eigen::Vector4d& currentQuaternionsToBaseFrame,
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
        // Initialize control vector to zero
        currentControlVector_ = Eigen::Vector3d::Zero( );
    }

    //! Destructor.
    ~ControlSystem( ) { }

    //! Attitude control system.
    void updateAttitudeController( const Eigen::Vector16d& currentEstimatedStateVector,
                                   const Eigen::Vector3d& currentMeasuredRotationalVelocityVector,
                                   const double navigationRefreshStepSize )
    {
        // Compute difference between current and commanded state
        Eigen::Vector4d currentErrorInEstimatedQuaternionState =
                computeErrorInEstimatedQuaternionState( currentEstimatedStateVector.segment( 6, 4 ),
                                                        computeCurrentCommandedQuaternionState( currentEstimatedStateVector ) );
        historyOfQuaternionStateErrors_.push_back( currentErrorInEstimatedQuaternionState.segment( 1, 3 ) );

//        std::cout << "Current error in quaternion: " << currentErrorInEstimatedQuaternionState.transpose( ) << std::endl;

        // Compute control vector based on control gains and error
        currentControlVector_ = proportionalGain_.cwiseProduct( currentErrorInEstimatedQuaternionState.segment( 1, 3 ) ) +
                integralGain_.cwiseProduct( numerical_quadrature::performExtendedSimpsonsQuadrature(
                                                navigationRefreshStepSize, historyOfQuaternionStateErrors_ ) ) +
                derivativeGain_.cwiseProduct( ( calculateQuaternionsDerivative(
                                                    currentEstimatedStateVector.segment( 6, 4 ),
                                                    currentMeasuredRotationalVelocityVector ) ).segment( 1, 3 ) );
        // only the imaginary part of the quaternion is used, since only three terms are needed to fully control the spacecraft
//        std::cout << "Current control vector: " << currentControlVector_.transpose( ) << std::endl;
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
    /*!
     *  Function to compute the current commanded quaternion to base frame. The body-fixed frame is assumed to correspond
     *  to the trajectory frame, where however, the y- and z-axes are inverted. Thus, the direction cosine matrix (DCM) describing
     *  the rotation from trajectory to inertial frame can be found first, then the second and third columns are inverted in
     *  sign, and finally the transpose is taken. The DCM is found by using the velocity and radial distance vector. The
     *  velocity (unit) vector corresponds directly to the x-axis of the trajectory frame, whereas the z-axis is computed by
     *  subtracting from the radial distance (unit) vector its projection on the x-axis. Then, the y-axis is determined via the
     *  right-hand rule (i.e., with the cross product). The commanded state thus corresponds to a state with zero angle of attack
     *  and angle of side-slip and with a bank angle of 180 degrees. Note that since the estimated state is used, the actual
     *  transformation can differ.
     *  \param currentEstimatedStateVector Current estimated state as provided by the navigation system.
     *  \return Quaternion representing the estimated inverse rotation from trajectory to inertial frame. Thus the commanded
     *      quaternion corresponds to a state with zero angle of attack and angle of side-slip and with a bank angle of 180 degrees.
     */
    Eigen::Vector4d computeCurrentCommandedQuaternionState( const Eigen::Vector16d& currentEstimatedStateVector )
    {
        // Declare direction cosine matrix
        Eigen::Matrix3d transformationFromTrajectoryToInertialFrame = Eigen::Matrix3d::Zero( );

        // Find the trajectory x-axis unit vector
        Eigen::Vector3d xUnitVector = currentEstimatedStateVector.segment( 3, 3 ).normalized( );
        transformationFromTrajectoryToInertialFrame.col( 0 ) = xUnitVector; // body-fixed (= trajectory)
//        std::cout << "x: " << xUnitVector.transpose( ) << std::endl;

        // Find trajectory z-axis unit vector
        Eigen::Vector3d zUnitVector = currentEstimatedStateVector.segment( 0, 3 ).normalized( );
        zUnitVector -= zUnitVector.dot( xUnitVector ) * xUnitVector;
        transformationFromTrajectoryToInertialFrame.col( 2 ) = - zUnitVector; // body-fixed (= -trajectory)
//        std::cout << "z: " << zUnitVector.transpose( ) << std::endl;

        // Find body-fixed y-axis unit vector
        transformationFromTrajectoryToInertialFrame.col( 1 ) = xUnitVector.cross( zUnitVector ); // body-fixed (= -trajectory)
//        std::cout << "DCM: " << std::endl << transformationFromTrajectoryToInertialFrame.transpose( ) << std::endl;

        // Transform DCM to quaternion and give output
        return linear_algebra::convertQuaternionToVectorFormat(
                    Eigen::Quaterniond( transformationFromTrajectoryToInertialFrame.transpose( ) ) );
        // the transpose is taken to return the inverse rotation from trajectory to inertial frame
    }

    //! Double denoting the proportional gain for the PID attitude controller.
    const Eigen::Vector3d proportionalGain_;

    //! Double denoting the integral gain for the PID attitude controller.
    const Eigen::Vector3d integralGain_;

    //! Double denoting the derivative gain for the PID attitude controller.
    const Eigen::Vector3d derivativeGain_;

    //! Vector denoting the current quaternion attitude correction.
    Eigen::Vector3d currentControlVector_;

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsoapsisManeuver_;

    //! History of errors in the estimated quaternion state.
    std::vector< Eigen::Vector3d > historyOfQuaternionStateErrors_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_CONTROL_SYSTEM_H
