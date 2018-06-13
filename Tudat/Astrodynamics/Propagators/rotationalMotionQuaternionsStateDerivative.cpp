/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/quaternionNormalizationMethod.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the matrix by which a quaternion vector is to be pre-multiplied to obtain this
//! quaternion's time-derivative.
Eigen::Matrix4d getQuaterionToQuaternionRateMatrix( const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
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

    return 0.5 * conversionMatrix;
}

//! Function to obtain the matrix by which an angular velocity vector is to be pre-multiplied to obtain the
//! quaternion's time-derivative.
Eigen::Matrix< double, 4, 3 > getAngularVelocityToQuaternionRateMatrix( const Eigen::Vector4d& quaternionVector )
{
    Eigen::Matrix< double, 4, 3 > conversionMatrix = Eigen::Matrix< double, 4, 3 >::Zero( );

    conversionMatrix( 0, 0 ) = -quaternionVector( 1 );
    conversionMatrix( 0, 1 ) = -quaternionVector( 2 );
    conversionMatrix( 0, 2 ) = -quaternionVector( 3 );

    conversionMatrix( 1, 0 ) = quaternionVector( 0 );
    conversionMatrix( 1, 1 ) = -quaternionVector( 3 );
    conversionMatrix( 1, 2 ) = quaternionVector( 2 );

    conversionMatrix( 2, 0 ) = quaternionVector( 3 );
    conversionMatrix( 2, 1 ) = quaternionVector( 0 );
    conversionMatrix( 2, 2 ) = -quaternionVector( 1 );

    conversionMatrix( 3, 0 ) = -quaternionVector( 2 );
    conversionMatrix( 3, 1 ) = quaternionVector( 1 );
    conversionMatrix( 3, 2 ) = quaternionVector( 0 );

    return 0.5 * conversionMatrix;
}

//! Function to obtain the time derivative of quaternions of body-fixed to inertial frame
Eigen::Vector4d calculateQuaternionsDerivative( const Eigen::Vector4d& currentQuaternionsToBaseFrame,
                                                const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
{
    switch ( QUATERNION_NORMALIZATION_METHOD )
    {
    case 0: // no internal normalization
        return getQuaterionToQuaternionRateMatrix( angularVelocityVectorInBodyFixedFrame ) * currentQuaternionsToBaseFrame;
    case 1: // normalization according to BOOK034 with rotational velocity
    {
        double gainFactor = angularVelocityVectorInBodyFixedFrame.norm( ) * ( 1.0 - currentQuaternionsToBaseFrame.norm( ) );
        return ( gainFactor * Eigen::Matrix4d::Identity( ) +
                 getQuaterionToQuaternionRateMatrix( angularVelocityVectorInBodyFixedFrame ) ) * currentQuaternionsToBaseFrame;
    }
    case 2: // normalization according to BOOK034 with constant factor
    {
        double gainFactor = 1e-6 * ( 1.0 - currentQuaternionsToBaseFrame.norm( ) );
        return ( gainFactor * Eigen::Matrix4d::Identity( ) +
                 getQuaterionToQuaternionRateMatrix( angularVelocityVectorInBodyFixedFrame ) ) * currentQuaternionsToBaseFrame;
    }
    case 3: // normalization of derivative accoding to ART073
    {
        Eigen::Vector4d quaternionDerivative = getQuaterionToQuaternionRateMatrix( angularVelocityVectorInBodyFixedFrame ) *
                currentQuaternionsToBaseFrame;
        return quaternionDerivative - currentQuaternionsToBaseFrame.dot( quaternionDerivative ) * currentQuaternionsToBaseFrame;
    }
    case 4: // normalization of quaternion and its derivative accoding to ART073
    {
        Eigen::Vector4d normalizedQuaternion = currentQuaternionsToBaseFrame.normalized( );
        Eigen::Vector4d quaternionDerivative = getQuaterionToQuaternionRateMatrix( angularVelocityVectorInBodyFixedFrame ) *
                normalizedQuaternion;
        return quaternionDerivative - normalizedQuaternion.dot( quaternionDerivative ) * normalizedQuaternion;
    }
    default:
        throw std::runtime_error( "Normalization method not recognized." );
    }
}

} // namespace propagators

} // namespace tudat
