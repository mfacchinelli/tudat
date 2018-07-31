#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

namespace tudat
{

namespace system_models
{

using namespace tudat::guidance_navigation_control;

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector12d onboardSystemModel( const double currentTime,
                                     const Eigen::Vector12d& currentEstimatedStateVector,
                                     const Eigen::Vector3d& currentEstimatedGravitationalTranslationalAccelerationVector,
                                     const Eigen::Vector3d& currentMeasuredNonGravitationalTranslationalAccelerationVector )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare state derivative vector
    Eigen::Vector12d currentStateDerivative = Eigen::Vector12d::Zero( );

    // Translational kinematics
    currentStateDerivative.segment( 0, 3 ) = currentEstimatedStateVector.segment( 3, 3 );

    // Translational dynamics
    Eigen::Vector3d currentActualNonGravitationalTranslationalAccelerationVector = removeErrorsFromInertialMeasurementUnitMeasurement(
                currentMeasuredNonGravitationalTranslationalAccelerationVector, currentEstimatedStateVector.segment( 6, 6 ) );
    currentStateDerivative.segment( 3, 3 ) = currentEstimatedGravitationalTranslationalAccelerationVector +
            currentActualNonGravitationalTranslationalAccelerationVector;

    // Give output
    return currentStateDerivative;
}

//! Function to model the onboard measurements based on the simplified onboard model.
Eigen::Vector1d onboardMeasurementModel( const double currentTime, const Eigen::Vector12d& currentEstimatedStateVector,
                                         const double planetaryRadius )
{
    TUDAT_UNUSED_PARAMETER( currentTime );

    // Declare output vector
    Eigen::Vector1d currentMeasurementVector;

    // Add translational acceleration
    currentMeasurementVector[ 0 ] = currentEstimatedStateVector.segment( 0, 3 ).norm( ) - planetaryRadius;

    // Return quaternion vector
    return currentMeasurementVector;
}

} // namespace thesis

} // namespace tudat
