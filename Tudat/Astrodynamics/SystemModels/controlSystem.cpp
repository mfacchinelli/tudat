#include "Tudat/Astrodynamics/SystemModels/controlSystem.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace system_models
{

//! Function to compute the error in the current quaternion state, based on the commanded quaternion.
Eigen::Vector4d computeErrorInEstimatedQuaternion( const Eigen::Vector4d& currentQuaternionToBaseFrame,
                                                   const Eigen::Vector4d& currentCommandedQuaternionToBaseFrame )
{
    // Define output vector
    Eigen::Vector4d currentErrorQuaternion = Eigen::Vector4d::Zero( );

    // Compute auxiliary matrix
    Eigen::Matrix< double, 3, 4 > auxiliaryMatrix = Eigen::Matrix< double, 3, 4 >::Zero( );
    auxiliaryMatrix.col( 0 ) = - currentCommandedQuaternionToBaseFrame.segment( 1, 3 );
    auxiliaryMatrix.rightCols( 3 ) = currentCommandedQuaternionToBaseFrame[ 0 ] * Eigen::Matrix3d::Identity( ) -
            linear_algebra::getCrossProductMatrix( currentCommandedQuaternionToBaseFrame.segment( 1, 3 ) );

    // Compute error quaternion
    currentErrorQuaternion.segment( 1, 3 ) = auxiliaryMatrix * currentQuaternionToBaseFrame;

    // Give output
    return currentErrorQuaternion;
}

//! Function to obtain the time derivative of quaternions of body-fixed to inertial frame
Eigen::Vector4d calculateQuaternionDerivative( const Eigen::Vector4d& currentQuaternionsToBaseFrame,
                                               const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
{
    Eigen::Matrix4d conversionMatrix = Eigen::Matrix4d::Zero( );
    conversionMatrix( 1, 0 ) = angularVelocityVectorInBodyFixedFrame( 0 );
    conversionMatrix( 2, 0 ) = angularVelocityVectorInBodyFixedFrame( 1 );
    conversionMatrix( 3, 0 ) = angularVelocityVectorInBodyFixedFrame( 2 );

    conversionMatrix( 2, 1 ) = -angularVelocityVectorInBodyFixedFrame( 2 );
    conversionMatrix( 3, 1 ) = angularVelocityVectorInBodyFixedFrame( 1 );

    conversionMatrix( 3, 2 ) = -angularVelocityVectorInBodyFixedFrame( 0 );

    conversionMatrix( 0, 1 ) = -angularVelocityVectorInBodyFixedFrame( 0 );
    conversionMatrix( 0, 2 ) = -angularVelocityVectorInBodyFixedFrame( 1 );
    conversionMatrix( 0, 3 ) = -angularVelocityVectorInBodyFixedFrame( 2 );

    conversionMatrix( 1, 2 ) = angularVelocityVectorInBodyFixedFrame( 2 );
    conversionMatrix( 1, 3 ) = -angularVelocityVectorInBodyFixedFrame( 1 );

    conversionMatrix( 2, 3 ) = angularVelocityVectorInBodyFixedFrame( 0 );

    conversionMatrix *= 0.5;

    return conversionMatrix * currentQuaternionsToBaseFrame;
}

//! Function to update the attitude controller.
void ControlSystem::updateAttitudeController( const Eigen::Vector4d& currentEstimatedQuaternion,
                                              const Eigen::Vector3d& currentMeasuredRotationalVelocityVector,
                                              const double navigationRefreshStepSize )
{
    // Determine commanded quaternion
    Eigen::Vector4d currentCommandedQuaternion = ( isSignToBeFlipped_ ? -1.0 : 1.0 ) * currentOrbitCommandedQuaternionState_;

    // Check if current quaternion signs match the previous time otherwise flip the sign
    bool doesSignHistoryMatch = currentCommandedQuaternion.isApprox( -previousCommandedQuaternion_, 1.0e-3 );
    if ( doesSignHistoryMatch )
    {
        currentCommandedQuaternion *= -1.0;
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
    // only the imaginary part of the quaternion is used, since only three terms are needed to fully control the spacecraft
//    std::cout << "Current control vector: " << currentControlVector_.transpose( ) << std::endl << std::endl;
    currentOrbitHistoryOfControlVectors_.push_back( currentControlVector_ );
}

//! Function to update the attitude controller.
/*!
 *  Function to update the attitude controller. The controller is based on the LQR (linear quadratic regulator).
 *  \param currentEstimatedAerodynamicAngles Current estimated aerodynamic angles.
 *  \param currentMeasuredRotationalVelocityVector Current measured rotational velocity vector.
 */
void ControlSystem::updateAttitudeControllerLQR( const Eigen::Vector3d& currentEstimatedAerodynamicAngles,
                                                 const Eigen::Vector3d& currentMeasuredRotationalVelocityVector )
{
    // Declare gain matrix
    Eigen::Matrix< double, 3, 6 > gainMatrix;
    gainMatrix << 2286.736461, 0.000000, 0.000000, 0.000000, 1860.175245, -0.000000,
            0.000000, -0.000000, 2639.026680, 0.000000, 1260.507149, 0.000000,
            0.000000, -0.000000, 630.253575, 420.169050, -0.000000, -0.000000;

    // Compute control vector
    Eigen::Vector6d inputVector;
    inputVector.segment( 0, 3 ) = currentMeasuredRotationalVelocityVector;
    inputVector.segment( 3, 3 ) = currentEstimatedAerodynamicAngles;
    currentControlVector_ = gainMatrix * inputVector / 1.0e5;
    std::cout << currentControlVector_.transpose( ) << std::endl;
}

} // namespace system_models

} // namespace tudat
