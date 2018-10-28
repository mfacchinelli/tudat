/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Facchinelli, M. (2018). Aerobraking Navigation, Guidance and Control.
 *          Master Thesis, Delft University of Technology.
 */

#ifndef TUDAT_CONTROL_SYSTEM_H
#define TUDAT_CONTROL_SYSTEM_H

#include <iostream>
#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Basics/utilityMacros.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace system_models
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
     *  \param proportionalGain Vector denoting the proportional gain for each of the imaginary quaternion elements.
     *  \param integralGain Vector denoting the integral gain for each of the imaginary quaternion elements.
     *  \param derivativeGain Vector denoting the derivative gain for each of the imaginary quaternion elements.
     */
    ControlSystem( const Eigen::Vector3d& proportionalGain,
                   const Eigen::Vector3d& integralGain,
                   const Eigen::Vector3d& derivativeGain ) :
        proportionalGain_( proportionalGain ), integralGain_( integralGain ), derivativeGain_( derivativeGain )
    {
        // Set values to their initial conditions
        currentControlVector_.setZero( );
        scheduledApsisManeuver_.setZero( );
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
     *  \param currentEstimatedQuaternion Current estimated quaternion to base frame.
     *  \param currentMeasuredRotationalVelocityVector Current measured rotational velocity vector.
     *  \param navigationRefreshStepSize Refresh step size of the navigation system.
     */
    void updateAttitudeController( const Eigen::Vector4d& currentEstimatedQuaternion,
                                   const Eigen::Vector3d& currentMeasuredRotationalVelocityVector,
                                   const double navigationRefreshStepSize )
    {
        // Determine commanded quaternion
        Eigen::Vector4d currentCommandedQuaternion = ( isSignToBeFlipped_ ? -1.0 : 1.0 ) * currentOrbitCommandedQuaternionState_;

        // Check if current quaternion signs match the previous time otherwise flip the sign
        bool doesSignHistoryMatch = currentCommandedQuaternion.isApprox( -previousCommandedQuaternion_, 1.0e-3 );
        if ( doesSignHistoryMatch )
        {
//            currentCommandedQuaternion *= -1.0;
        }
        isSignToBeFlipped_ = isSignToBeFlipped_ || doesSignHistoryMatch;
        previousCommandedQuaternion_ = currentCommandedQuaternion;

        // Compute difference between current and commanded state
        Eigen::Vector4d currentErrorInEstimatedQuaternion =
                computeErrorInEstimatedQuaternion( currentEstimatedQuaternion, currentCommandedQuaternion );
        historyOfQuaternionErrors_.push_back( currentErrorInEstimatedQuaternion.segment( 1, 3 ) );

        // Compute difference between current and commanded derivative
        Eigen::Vector4d currentErrorInEstimatedQuaternionDerivative =
                computeErrorInEstimatedQuaternion( calculateQuaternionDerivative( currentEstimatedQuaternion,
                                                                                  currentMeasuredRotationalVelocityVector ),
                                                   currentCommandedQuaternion );

//        std::cout << "Current estimated quaternion: " <<
//                     currentEstimatedQuaternion.transpose( ) << std::endl
//                  << "Current estimated derivative: "  <<
//                     calculateQuaternionDerivative( currentEstimatedQuaternion,
//                                                    currentMeasuredRotationalVelocityVector ).transpose( ) << std::endl
//                  << "Commanded state: " << currentCommandedQuaternion.transpose( ) << std::endl
//                  << "Proportional: " <<
//                     proportionalGain_.cwiseProduct( currentErrorInEstimatedQuaternion.segment( 1, 3 ) ).transpose( ) << std::endl
//                  << "Integral: " <<
//                     integralGain_.cwiseProduct( numerical_quadrature::performExtendedSimpsonsQuadrature(
//                                                     navigationRefreshStepSize, historyOfQuaternionErrors_ ) ).transpose( ) << std::endl
//                  << "Derivative: " <<
//                     derivativeGain_.cwiseProduct( currentErrorInEstimatedQuaternionDerivative.segment( 1, 3 ) ).transpose( ) << std::endl;

        // Compute control vector based on control gains and error
        currentControlVector_ = - ( proportionalGain_.cwiseProduct( currentErrorInEstimatedQuaternion.segment( 1, 3 ) ) +
                                    integralGain_.cwiseProduct( numerical_quadrature::performExtendedSimpsonsQuadrature(
                                                                    navigationRefreshStepSize, historyOfQuaternionErrors_ ) ) +
                                    derivativeGain_.cwiseProduct( currentErrorInEstimatedQuaternionDerivative.segment( 1, 3 ) ) );
        currentOrbitHistoryOfControlVectors_.push_back( currentControlVector_ );
        // only the imaginary part of the quaternion is used, since only three terms are needed to fully control the spacecraft
//        std::cout << "Current control vector: " << currentControlVector_.transpose( ) << std::endl << std::endl;
    }

    //! Function to update the orbit controller with the scheduled apoapsis maneuver, computed by the guidance system.
    void updateOrbitController( const Eigen::Vector3d& scheduledApsisManeuver )
    {
        scheduledApsisManeuver_ = scheduledApsisManeuver;
    }

    //! Function to retireve current control vector for attitude.
    Eigen::Vector3d getCurrentAttitudeControlVector( ) { return currentControlVector_; }

    //! Function to retirieve the apoapsis maneuver.
    Eigen::Vector3d getScheduledApsisManeuver( ) { return scheduledApsisManeuver_; }

    //! Function to retrieve history of commanded quaternions.
    std::map< unsigned int, Eigen::Vector4d > getHistoryOfCommandedQuaternions( ) { return historyOfCommandedQuaternionStates_; }

    //! Function to retrieve history of control vectors for the current orbit.
    std::vector< Eigen::Vector3d > getCurrentOrbitHistoryOfControlVectors( ) { return currentOrbitHistoryOfControlVectors_; }

    //! Function to set the value of the current orbit counter.
    void setCurrentOrbitCounter( const unsigned int currentOrbitCounter )
    {
        currentOrbitCounter_ = currentOrbitCounter;
    }

    //! Function to set the commanded quaternion based on periapsis conditions.
    void setCommandedQuaternionBasedOnPeriapsisConditions( const Eigen::Vector6d& estimatedTranslationalStateAtPeriapsis )
    {
        currentOrbitCommandedQuaternionState_ = computeCurrentCommandedQuaternion( estimatedTranslationalStateAtPeriapsis );
    }

    //! Clear history of control vectors for the current orbit.
    void clearCurrentOrbitControlHistory( )
    {
        currentOrbitHistoryOfControlVectors_.clear( );
    }

    //! Function to revert to the previous time step.
    /*!
     *  Function to revert to the previous time step. This function is run if the current propagation needs to be stopped, since
     *  the current time will be run the next time the GNC system is called.
     *  \param currentTime Double denoting the current time, i.e., the instant that has to be discarded.
     */
    void revertToPreviousTimeStep( const double currentTime )
    {
        TUDAT_UNUSED_PARAMETER( currentTime );
    }

private:

    //! Function to compute the current commanded quaternion to base frame.
    /*!
     *  Function to compute the current commanded quaternion to base frame. The body-fixed frame is assumed to correspond
     *  to the trajectory frame. Thus, the direction cosine matrix (DCM) describing the rotation from trajectory to inertial
     *  frame can be found and taken as full rotation from body-fixed to inertial. The DCM is found by using the velocity and
     *  radial distance vector. The airspeed (unit) vector corresponds directly to the x-axis of the trajectory frame, whereas
     *  the z-axis is computed by subtracting from the radial distance (unit) vector its projection on the x-axis. Then, the
     *  y-axis is determined via the right-hand rule (i.e., with the cross product). The commanded state thus corresponds to a
     *  state with zero angle of attack, angle of side-slip and bank angle. Note that since the estimated state is used, the actual
     *  transformation can differ.
     *  \param estimatedTranslationalStateAtPeriapsis Estimated translational Cartesian state at periapsis.
     *  \return Quaternion representing the estimated rotation from trajectory to inertial frame. Thus the commanded
     *      quaternion corresponds to a state with zero angle of attack, angle of side-slip and bank angle.
     */
    Eigen::Vector4d computeCurrentCommandedQuaternion( const Eigen::Vector6d& estimatedTranslationalStateAtPeriapsis )
    {
        // Declare eventual ouptut vector
        Eigen::Vector4d commandedQuaternion;

        // Declare direction cosine matrix
        Eigen::Matrix3d transformationFromInertialToTrajectoryFrame;

        // Find the trajectory x-axis unit vector
        Eigen::Vector3d currentRadialVector = estimatedTranslationalStateAtPeriapsis.segment( 0, 3 );
        Eigen::Vector3d xUnitVector = ( estimatedTranslationalStateAtPeriapsis.segment( 3, 3 ) -
                                        7.07763225880808e-05 * Eigen::Vector3d::UnitZ( ).cross( currentRadialVector ) ).normalized( );
        transformationFromInertialToTrajectoryFrame.col( 0 ) = xUnitVector;

        // Find trajectory z-axis unit vector
        Eigen::Vector3d zUnitVector = currentRadialVector.normalized( );
        zUnitVector -= zUnitVector.dot( xUnitVector ) * xUnitVector;
        transformationFromInertialToTrajectoryFrame.col( 2 ) = zUnitVector;

        // Find body-fixed y-axis unit vector
        transformationFromInertialToTrajectoryFrame.col( 1 ) = zUnitVector.cross( xUnitVector );

        // Transform DCM to quaternion and give output
        commandedQuaternion = linear_algebra::convertQuaternionToVectorFormat(
                    Eigen::Quaterniond( transformationFromInertialToTrajectoryFrame ) );
        historyOfCommandedQuaternionStates_[ currentOrbitCounter_ ] = commandedQuaternion;
        return commandedQuaternion;
        // note that due to the different rotation convention in Eigen, transforming the DCM with the constructor above
        // automatically takes the inverse of the rotation
    }

    //! Double denoting the proportional gain for the PID attitude controller.
    const Eigen::Vector3d proportionalGain_;

    //! Double denoting the integral gain for the PID attitude controller.
    const Eigen::Vector3d integralGain_;

    //! Double denoting the derivative gain for the PID attitude controller.
    const Eigen::Vector3d derivativeGain_;

    //! Vector denoting the previous signs of quaternions.
    Eigen::Vector4d previousCommandedQuaternion_;

    //! Boolean denoting whether the sign of the quaternion needs to be flipped.
    bool isSignToBeFlipped_;

    //! Integer denoting the value of the current orbit.
    unsigned int currentOrbitCounter_;

    //! Vector denoting the commanded quaternion based on periapsis conditions.
    Eigen::Vector4d currentOrbitCommandedQuaternionState_;

    //! Vector denoting the current quaternion attitude correction.
    Eigen::Vector3d currentControlVector_;

    //! Vector denoting the velocity change scheduled to be applied at apo- or periapsis.
    Eigen::Vector3d scheduledApsisManeuver_;

    //! History of commanded quaternion states.
    std::map< unsigned int, Eigen::Vector4d > historyOfCommandedQuaternionStates_;

    //! History of errors in the estimated quaternion state.
    std::vector< Eigen::Vector3d > historyOfQuaternionErrors_;

    //! History of control torque vectors for current orbit.
    std::vector< Eigen::Vector3d > currentOrbitHistoryOfControlVectors_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_CONTROL_SYSTEM_H
