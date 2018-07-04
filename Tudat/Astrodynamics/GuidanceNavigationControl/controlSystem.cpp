#include "Tudat/Astrodynamics/GuidanceNavigationControl/controlSystem.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace guidance_navigation_control
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

} // namespace navigation

} // namespace tudat
