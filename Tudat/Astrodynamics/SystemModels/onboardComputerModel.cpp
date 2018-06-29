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
                                     const Eigen::Vector3d& currentRotationalVelocityVector )
{
    // Declare state derivative vector
    Eigen::Vector16d currentStateDerivative = Eigen::Vector16d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentStateVector.segment( 3, 3 );

    // Translational dynamics
    currentStateDerivative.segment( 3, 3 ) = currentTranslationalAccelerationVector;

    // Rotational kinematics
    Eigen::Vector3d uncorruptedCurrentRotationalVelocityVector =
            ( Eigen::Matrix3d::Identity( ) - currentStateVector.segment( 13, 3 ).asDiagonal( ) ) *
            ( currentRotationalVelocityVector - currentStateVector.segment( 10, 3 ) );
    currentStateDerivative.segment( 6, 4 ) = propagators::calculateQuaternionsDerivative( currentStateVector.segment( 6, 4 ),
                                                                                          uncorruptedCurrentRotationalVelocityVector );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector4d onboardMeasurementModel( const double currentTime, const Eigen::Vector16d& currentStateVector )
{
    // Return quaternion vector
    return currentStateVector.segment( 6, 4 );
}

} // namespace thesis

} // namespace tudat
