#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

namespace tudat
{

namespace system_models
{

using namespace guidance_navigation_control;

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector12d onboardSystemModel(
        const double currentTime, const Eigen::Vector12d& currentEstimatedStateVector,
        const boost::function< Eigen::Vector3d( const Eigen::Vector6d& ) >& currentEstimatedTranslationalAccelerationVectorFunction )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare state derivative vector
    Eigen::Vector12d currentStateDerivative = Eigen::Vector12d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentEstimatedStateVector.segment( 3, 3 );

    // Translational dynamics
    currentStateDerivative.segment( 3, 3 ) =
            currentEstimatedTranslationalAccelerationVectorFunction( currentEstimatedStateVector.segment( 0, 6 ) );

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector3d onboardMeasurementModel(
        const double currentTime, const Eigen::Vector12d& currentEstimatedStateVector,
        const boost::function< Eigen::Vector3d( const Eigen::Vector6d& ) >& currenstEstimatedNonGravitationalAccelerationFunction )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::Vector3d currentMeasurementVector;

    // Add translational acceleration
    currentMeasurementVector = currentEstimatedStateVector.segment( 6, 3 ) +
            ( Eigen::Matrix3d::Identity( ) + Eigen::Matrix3d( currentEstimatedStateVector.segment( 9, 3 ).asDiagonal( ) ) ) *
            currenstEstimatedNonGravitationalAccelerationFunction( currentEstimatedStateVector.segment( 0, 6 ) );

    // Return quaternion vector
    return currentMeasurementVector;
}

} // namespace thesis

} // namespace tudat
