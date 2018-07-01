#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

namespace tudat
{

namespace system_models
{

using namespace tudat::guidance_navigation_control;

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector16d onboardSystemModel( const double currentTime,
                                     const Eigen::Vector16d& currentStateVector,
                                     const Eigen::VectorXd& currentControlVector,
                                     const Eigen::Vector3d& currentTranslationalAccelerationVector,
                                     const Eigen::Vector3d& currentMeasuredRotationalVelocityVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare state derivative vector
    Eigen::Vector16d currentStateDerivative = Eigen::Vector16d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentStateVector.segment( 3, 3 );

    // Translational dynamics
    currentStateDerivative.segment( 3, 3 ) = currentTranslationalAccelerationVector;

    // Rotational kinematics
    Eigen::Vector3d uncorruptedCurrentRotationalVelocityVector =
            ( Eigen::Matrix3d::Identity( ) - currentStateVector.segment( 13, 3 ).asDiagonal( ) ) * // binomial approximation
            ( currentMeasuredRotationalVelocityVector - currentStateVector.segment( 10, 3 ) ); // + currentControlVector
    currentStateDerivative.segment( 6, 4 ) = propagators::calculateQuaternionsDerivative( currentStateVector.segment( 6, 4 ),
                                                                                          uncorruptedCurrentRotationalVelocityVector );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector7d onboardMeasurementModel( const double currentTime, const Eigen::Vector16d& currentStateVector,
                                         const Eigen::Vector3d& currentTranslationalAccelerationVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::Vector7d currentMeasurementVector;

    // Add translational acceleration
    currentMeasurementVector.segment( 0, 3 ) = currentTranslationalAccelerationVector;

    // Add rotational attitude
    currentMeasurementVector.segment( 3, 4 ) = currentStateVector.segment( 6, 4 );

    // Return quaternion vector
    return currentMeasurementVector;
}

} // namespace thesis

} // namespace tudat
