#include "Tudat/Astrodynamics/GuidanceNavigationControl/controlSystem.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to compute the error in the current quaternion state, based on the commanded quaternion.
Eigen::Vector4d computeErrorInEstimatedQuaternionState( const Eigen::Vector4d& currentQuaternionToBaseFrame,
                                                        const Eigen::Vector4d& currentCommandedQuaternionToBaseFrame )
{
    // Define output vector
    Eigen::Vector4d currentErrorQuaternion = Eigen::Vector4d::Zero( );

    // Compute auxiliary matrix
    Eigen::Matrix< double, 3, 4 > auxiliaryMatrix = Eigen::Matrix< double, 3, 4 >::Zero( );
    auxiliaryMatrix.block( 0, 0, 3, 1 ) = - currentCommandedQuaternionToBaseFrame.segment( 1, 3 );
    auxiliaryMatrix.block( 0, 1, 3, 3 ) = currentCommandedQuaternionToBaseFrame[ 0 ] * Eigen::Matrix3d::Identity( ) -
            linear_algebra::getCrossProductMatrix( currentCommandedQuaternionToBaseFrame.segment( 1, 3 ) );

    // Compute error quaternion
    currentErrorQuaternion.segment( 1, 3 ) = auxiliaryMatrix * currentQuaternionToBaseFrame;
}

} // namespace navigation

} // namespace tudat
