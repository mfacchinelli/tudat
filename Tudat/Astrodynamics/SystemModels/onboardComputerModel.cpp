#include "Tudat/Astrodynamics/SystemModels/onboardComputerModel.h"

#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"

namespace tudat
{

namespace system_models
{

using namespace tudat::guidance_navigation_control;

//! Function to remove the error in gyroscope measurement based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromGyroscopeMeasurement( const Eigen::Vector3d& currentGyroscopeMeasurement,
                                                      const Eigen::Vector16d& currentEstimatedStateVector )
{
    return ( Eigen::Matrix3d::Identity( ) -
             Eigen::Matrix3d( currentEstimatedStateVector.segment( 13, 3 ).asDiagonal( ) ) ) * // binomial approximation
           ( currentGyroscopeMeasurement - currentEstimatedStateVector.segment( 10, 3 ) );
}

//! Function to model the onboard system dynamics based on the simplified onboard model.
Eigen::Vector16d onboardSystemModel( const double currentTime,
                                     const Eigen::Vector16d& currentEstimatedStateVector,
                                     const Eigen::Vector3d& currentControlVector,
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
    std::cout << "Full: " << currentEstimatedTranslationalAccelerationVector.transpose( ) << std::endl;

    // Rotational kinematics
    Eigen::Vector3d actualCurrentRotationalVelocityVector = removeErrorsFromGyroscopeMeasurement(
                currentMeasuredRotationalVelocityVector, currentEstimatedStateVector ) + currentControlVector;
    currentStateDerivative.segment( 6, 4 ) = propagators::calculateQuaternionsDerivative( currentEstimatedStateVector.segment( 6, 4 ),
                                                                                          actualCurrentRotationalVelocityVector );

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
    std::cout << "Non-grav.: " << currentEstimatedNonGravitationalTranslationalAccelerationVector.transpose( ) << std::endl;

    // Add rotational attitude
    currentMeasurementVector.segment( 3, 4 ) = currentEstimatedStateVector.segment( 6, 4 );

    // Return quaternion vector
    return currentMeasurementVector;
}

} // namespace thesis

} // namespace tudat
