#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"

namespace tudat
{

namespace system_models
{

//! Function to compute the scale-misalignment matrix from the scale and misalignment vectors.
Eigen::Matrix3d computeScaleMisalignmentMatrix( const Eigen::Vector3d& scaleFactorVector,
                                                const Eigen::Vector6d& misalignmentVector )
{
    // Compute misalignment matrix
    Eigen::Matrix3d misalignmentMatrix = Eigen::Matrix3d::Zero( );
    misalignmentMatrix( 0, 1 ) = misalignmentVector[ 0 ];
    misalignmentMatrix( 0, 2 ) = misalignmentVector[ 1 ];
    misalignmentMatrix( 1, 0 ) = misalignmentVector[ 2 ];
    misalignmentMatrix( 1, 2 ) = misalignmentVector[ 3 ];
    misalignmentMatrix( 2, 0 ) = misalignmentVector[ 4 ];
    misalignmentMatrix( 2, 1 ) = misalignmentVector[ 5 ];

    // Give output
    return Eigen::Matrix3d::Identity( ) + scaleFactorVector.asDiagonal( ) + misalignmentMatrix;
}

//! Function to add an inertial measurement unit to the spacecraft set of instruments.
void NavigationInstrumentsModel::addInertialMeasurementUnit( const Eigen::Vector3d& accelerometerBias,
                                                             const Eigen::Vector3d& accelerometerScaleFactor,
                                                             const Eigen::Vector6d& accelerometerMisalignment,
                                                             const Eigen::Vector3d& accelerometerAccuracy,
                                                             const Eigen::Vector3d& gyroscopeBias,
                                                             const Eigen::Vector3d& gyroscopeScaleFactor,
                                                             const Eigen::Vector6d& gyroscopeMisalignment,
                                                             const Eigen::Vector3d& gyroscopeAccuracy )
{
    // Check whether an inertial measurement unit is already present
    if ( !inertialMeasurementUnitAdded_ )
    {
        // Inertial measurement unit has been created
        inertialMeasurementUnitAdded_ = true;

        // Generate random noise distribution
        generateInertialMeasurementUnitRandomNoiseDistribution( accelerometerAccuracy, gyroscopeAccuracy );

        // Create function for computing corrupted translational accelerations
        inertialMeasurementUnitTranslationalAccelerationFunction_ = boost::bind(
                    &NavigationInstrumentsModel::getCurrentTranslationalAcceleration, this, accelerometerBias,
                    computeScaleMisalignmentMatrix( accelerometerScaleFactor, accelerometerMisalignment ) );

        // Create function for computing corrupted rotational velocity
        inertialMeasurementUnitRotationalVelocityFunction_ = boost::bind(
                    &NavigationInstrumentsModel::getCurrentRotationalVelocity, this, gyroscopeBias,
                    computeScaleMisalignmentMatrix( gyroscopeScaleFactor, gyroscopeMisalignment ) );
    }
    else
    {
        throw std::runtime_error( "Error in creation of inertial measurement unit for body " + spacecraftName_ +
                                  ". An IMU is already present." );
    }
}

//! Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
void NavigationInstrumentsModel::addStarTracker( const unsigned int numberOfStarTrackers,
                                                 const Eigen::Vector3d& starTrackerAccuracy )
{
    // Check whether a star tracker system is already present
    if ( !starTrackerAdded_ )
    {
        if ( numberOfStarTrackers == 2 )
        {
            // Star tracker has been created
            starTrackerAdded_ = true;

            // Generate random noise distribution
            generateStarTrackerRandomNoiseDistribution( starTrackerAccuracy );

            // Create function for computing corrupted spacecraft orientation
            starTrackerOrientationFunction_ = boost::bind(
                        &NavigationInstrumentsModel::getCurrentAttitude, this );
        }
        else
        {
            throw std::runtime_error( "Error in creation of star tracker system for body " + spacecraftName_ +
                                      ". Only a system of two orthogonal star trackers is supported." );
        }
    }
    else
    {
        throw std::runtime_error( "Error in creation of star tracker system for body " + spacecraftName_ +
                                  ". A star tracker system is already present." );
    }
}

//! Function to generate the noise distributions for the inertial measurement unit.
void NavigationInstrumentsModel::generateInertialMeasurementUnitRandomNoiseDistribution(
        const Eigen::Vector3d& accelerometerAccuracy,
        const Eigen::Vector3d& gyroscopeAccuracy )
{
    using namespace tudat::statistics;

    // Create accelerometer noise distribution
    for ( unsigned int i = 0; i < 3; i++ )
    {
        if ( accelerometerAccuracy[ i ] != 0.0 )
        {
            accelerometerNoiseDistribution_.push_back(
                        createBoostContinuousRandomVariableGenerator(
                            normal_boost_distribution, { 0.0, accelerometerAccuracy[ i ] / 3.0 }, i ) );
        }
        else
        {
            accelerometerNoiseDistribution_.push_back( NULL );
        }
    }

    // Create gyroscope noise distribution
    for ( unsigned int i = 0; i < 3; i++ )
    {
        if ( gyroscopeAccuracy[ i ] != 0.0 )
        {
            gyroscopeNoiseDistribution_.push_back(
                        createBoostContinuousRandomVariableGenerator(
                            normal_boost_distribution, { 0.0, gyroscopeAccuracy[ i ] / 3.0 },
                            accelerometerAccuracy.rows( ) + i ) );
        }
        else
        {
            gyroscopeNoiseDistribution_.push_back( NULL );
        }
    }
}

//! Function to generate the noise distributions for the star trackers.
void NavigationInstrumentsModel::generateStarTrackerRandomNoiseDistribution(
        const Eigen::Vector4d& starTrackerAccuracy )
{
    using namespace tudat::statistics;

    // Create accelerometer noise distribution
    for ( unsigned int i = 0; i < 4; i++ )
    {
        if ( starTrackerAccuracy[ i ] != 0.0 )
        {
            starTrackerNoiseDistribution_.push_back(
                        createBoostContinuousRandomVariableGenerator(
                            normal_boost_distribution, { 0.0, starTrackerAccuracy[ i ] / 3.0 }, i ) );
        }
        else
        {
            starTrackerNoiseDistribution_.push_back( NULL );
        }
    }
}

} // namespace system_models

} // namespace tudat