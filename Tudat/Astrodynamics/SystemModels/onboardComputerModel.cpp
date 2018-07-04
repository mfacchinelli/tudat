#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

namespace tudat
{

namespace system_models
{

using namespace tudat::guidance_navigation_control;

//! Function to remove the error in gyroscope measurement based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromInertialMeasurementUnitMeasurement( const Eigen::Vector3d& currentInertialMeasurementUnitMeasurement,
                                                                    const Eigen::Vector16d& currentEstimatedStateVector )
{
    return ( Eigen::Matrix3d::Identity( ) -
             Eigen::Matrix3d( currentEstimatedStateVector.segment( 13, 3 ).asDiagonal( ) ) ) * // binomial approximation
            ( currentInertialMeasurementUnitMeasurement - currentEstimatedStateVector.segment( 10, 3 ) );
}

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector16d onboardSystemModel( const double currentTime,
                                     const Eigen::Vector16d& currentEstimatedStateVector,
                                     const Eigen::Vector3d& currentControlVector,
                                     const Eigen::Vector3d& currentEstimatedTranslationalAccelerationVector,
                                     const Eigen::Vector3d& currentMeasuredRotationalVelocityVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );
    TUDAT_UNUSED_PARAMETER( currentControlVector );

    // Declare state derivative vector
    Eigen::Vector16d currentStateDerivative = Eigen::Vector16d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentEstimatedStateVector.segment( 3, 3 );

    // Translational dynamics
    currentStateDerivative.segment( 3, 3 ) = currentEstimatedTranslationalAccelerationVector;

    // Invert quaternion to get rotation to global (base) frame
    Eigen::Vector4d currentQuaternionToGlobalFrame = currentEstimatedStateVector.segment( 6, 4 );
    currentQuaternionToGlobalFrame.segment( 1, 3 ) *= -1.0;

    // Rotational kinematics
    Eigen::Vector3d currentActualRotationalVelocityVector = removeErrorsFromInertialMeasurementUnitMeasurement(
                currentMeasuredRotationalVelocityVector, currentEstimatedStateVector );
    currentStateDerivative.segment( 6, 4 ) = propagators::calculateQuaternionDerivative(
                currentQuaternionToGlobalFrame, currentActualRotationalVelocityVector );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector7d onboardMeasurementModel( const double currentTime, const Eigen::Vector16d& currentEstimatedStateVector,
                                         const Eigen::Vector3d& currentEstimatedNonGravitationalTranslationalAccelerationVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::Vector7d currentMeasurementVector;

    // Add translational acceleration
    currentMeasurementVector.segment( 0, 3 ) = currentEstimatedNonGravitationalTranslationalAccelerationVector;

    // Add rotational attitude
    currentMeasurementVector.segment( 3, 4 ) = currentEstimatedStateVector.segment( 6, 4 );

    // Return quaternion vector
    return currentMeasurementVector;
}

} // namespace thesis

} // namespace tudat
