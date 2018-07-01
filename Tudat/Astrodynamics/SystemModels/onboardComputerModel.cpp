#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

namespace tudat
{

namespace system_models
{

using namespace tudat::guidance_navigation_control;

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector16d onboardSystemModel( const double currentTime,
                                     const Eigen::Vector16d& currentEstimatedStateVector,
                                     const Eigen::VectorXd& currentControlVector,
                                     const Eigen::Vector3d& currentEstimatedTranslationalAccelerationVector,
                                     const Eigen::Vector3d& currentMeasuredRotationalVelocityVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare state derivative vector
    Eigen::Vector16d currentStateDerivative = Eigen::Vector16d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentEstimatedStateVector.segment( 3, 3 );

    // Translational dynamics
    currentStateDerivative.segment( 3, 3 ) = currentEstimatedTranslationalAccelerationVector;

    // Rotational kinematics
    Eigen::Vector3d uncorruptedCurrentRotationalVelocityVector =
            ( Eigen::Matrix3d::Identity( ) -
              Eigen::Matrix3d( currentEstimatedStateVector.segment( 13, 3 ).asDiagonal( ) ) ) * // binomial approximation
            ( currentMeasuredRotationalVelocityVector - currentEstimatedStateVector.segment( 10, 3 ) ); // + currentControlVector
    currentStateDerivative.segment( 6, 4 ) = propagators::calculateQuaternionsDerivative( currentEstimatedStateVector.segment( 6, 4 ),
                                                                                          uncorruptedCurrentRotationalVelocityVector );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector7d onboardMeasurementModel( const double currentTime, const Eigen::Vector16d& currentEstimatedStateVector,
                                         const Eigen::Vector3d& currentEstimatedTranslationalAccelerationVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::Vector7d currentMeasurementVector;

    // Add translational acceleration
    currentMeasurementVector.segment( 0, 3 ) = currentEstimatedTranslationalAccelerationVector;

    // Add rotational attitude
    currentMeasurementVector.segment( 3, 4 ) = currentEstimatedStateVector.segment( 6, 4 );

    // Return quaternion vector
    return currentMeasurementVector;
}

} // namespace thesis

} // namespace tudat
