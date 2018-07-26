#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

namespace tudat
{

namespace system_models
{

using namespace tudat::guidance_navigation_control;

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector22d onboardSystemModel( const double currentTime,
                                     const Eigen::Vector22d& currentEstimatedStateVector,
                                     const Eigen::Vector3d& currentControlVector,
                                     const Eigen::Vector3d& currentEstimatedGravitationalTranslationalAccelerationVector,
                                     const Eigen::Vector3d& currentMeasuredNonGravitationalTranslationalAccelerationVector,
                                     const Eigen::Vector3d& currentMeasuredRotationalVelocityVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );
    TUDAT_UNUSED_PARAMETER( currentControlVector );

    // Declare state derivative vector
    Eigen::Vector22d currentStateDerivative = Eigen::Vector22d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentEstimatedStateVector.segment( 3, 3 );

    // Translational dynamics
    Eigen::Vector3d currentActualNonGravitationalTranslationalAccelerationVector = removeErrorsFromInertialMeasurementUnitMeasurement(
                currentMeasuredNonGravitationalTranslationalAccelerationVector, currentEstimatedStateVector.segment( 10, 6 ) );
    currentStateDerivative.segment( 3, 3 ) = ( linear_algebra::convertVectorToQuaternionFormat(
                                                   currentEstimatedStateVector.segment( 6, 4 ) ).toRotationMatrix( ).transpose( ) *
                                               currentActualNonGravitationalTranslationalAccelerationVector ) +
            currentEstimatedGravitationalTranslationalAccelerationVector;
    // transpose is taken due to the different definition of quaternion in Eigen

    // Rotational kinematics
    Eigen::Vector3d currentActualRotationalVelocityVector = removeErrorsFromInertialMeasurementUnitMeasurement(
                currentMeasuredRotationalVelocityVector, currentEstimatedStateVector.segment( 16, 6 ) );
    currentStateDerivative.segment( 6, 4 ) = propagators::calculateQuaternionDerivative(
                currentEstimatedStateVector.segment( 6, 4 ).normalized( ), currentActualRotationalVelocityVector );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector5d onboardMeasurementModel( const double currentTime, const Eigen::Vector22d& currentEstimatedStateVector,
                                         const double planetaryRadius )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::Vector5d currentMeasurementVector;

    // Add translational acceleration
    currentMeasurementVector[ 0 ] = currentEstimatedStateVector.segment( 0, 3 ).norm( ) - planetaryRadius;

    // Add rotational attitude
    currentMeasurementVector.segment( 1, 4 ) = currentEstimatedStateVector.segment( 6, 4 ).normalized( );

    // Return quaternion vector
    return currentMeasurementVector;
}

} // namespace thesis

} // namespace tudat
