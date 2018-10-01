#include "Tudat/Astrodynamics/SystemModels/instrumentsModel.h"

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
    return Eigen::Matrix3d::Identity( ) + Eigen::Matrix3d( scaleFactorVector.asDiagonal( ) ) + misalignmentMatrix;
}

//! Function to add a small angle uncertainty to a quaternion vector.
Eigen::Vector4d sumQuaternionUncertainty( const Eigen::Vector4d& quaternionVector,
                                          const Eigen::Vector3d& equivalentVectorUncertainty )
{
    // Find quaternion representing the uncertainty in the angle vector
    // Here the small angle approximation is used, such that the transformation can be simplified
    Eigen::Vector4d quaternionUncertainty;
    quaternionUncertainty[ 0 ] = 1.0;
    quaternionUncertainty.segment( 1, 3 ) = 0.5 * equivalentVectorUncertainty;

    // Add the two vectors with the conventional quaternion product function
    return linear_algebra::quaternionProduct( quaternionVector, quaternionUncertainty );
}

//! Function to add an inertial measurement unit to the spacecraft set of instruments.
void InstrumentsModel::addInertialMeasurementUnit( const Eigen::Vector3d& accelerometerBias,
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
        std::cout << "Acc. Bias: " << accelerometerBias.transpose( ) << ". Scale: " << accelerometerScaleFactor.transpose( ) << std::endl;
        std::cout << "Gyr. Bias: " << gyroscopeBias.transpose( ) << ". Scale: " << gyroscopeScaleFactor.transpose( ) << std::endl;

        // Create function for computing corrupted translational accelerations
        inertialMeasurementUnitTranslationalAccelerationFunction_ = boost::bind(
                    &InstrumentsModel::getCurrentTranslationalAcceleration, this, accelerometerBias,
                    computeScaleMisalignmentMatrix( accelerometerScaleFactor, accelerometerMisalignment ) );

        // Create function for computing corrupted rotational velocity
        inertialMeasurementUnitRotationalVelocityFunction_ = boost::bind(
                    &InstrumentsModel::getCurrentRotationalVelocity, this, gyroscopeBias,
                    computeScaleMisalignmentMatrix( gyroscopeScaleFactor, gyroscopeMisalignment ) );
    }
    else
    {
        throw std::runtime_error( "Error in creation of inertial measurement unit for body " + spacecraftName_ +
                                  ". An IMU is already present." );
    }
}

//! Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
void InstrumentsModel::addStarTracker( const unsigned int numberOfStarTrackers,
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
            starTrackerOrientationFunction_ = boost::bind( &InstrumentsModel::getCurrentAttitude, this );
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

//! Function to add an altimeter to the spacecraft set of instruments.
void InstrumentsModel::addAltimeter( const std::vector< Eigen::Vector3d >& altimeterPointingDirectionInAltimeterFrame,
                                     const reference_frames::AerodynamicsReferenceFrames altimeterFrame,
                                     const std::pair< double, double >& altitudeRange,
                                     const boost::function< double( double ) >& altimeterAccuracyAsAFunctionOfAltitude )
{
    // Check whether an altimeter is already present
    if ( !altimeterAdded_ )
    {
        // Altimeter has been created
        altimeterAdded_ = true;

        // Generate random noise distribution
        generateAltimeterRandomNoiseDistribution( 1.0 );

        // Resize altimeter output
        currentAltitude_.resize( altimeterPointingDirectionInAltimeterFrame.size( ) );

        // Create function for computing corrupted spacecraft orientation
        altimeterFunction_ = boost::bind( &InstrumentsModel::getCurrentAltitude, this, altimeterPointingDirectionInAltimeterFrame,
                                          altimeterFrame, altitudeRange, altimeterAccuracyAsAFunctionOfAltitude );
    }
    else
    {
        throw std::runtime_error( "Error in creation of altimeter for body " + spacecraftName_ +
                                  ". An altimeter is already present." );
    }
}

//! Function to add a Deep Space Network measurement system.
void InstrumentsModel::addDeepSpaceNetwork( const double positionAccuracy,
                                            const double velocityAccuracy,
                                            const double lightTimeAccuracy )
{
    // Check whether an altimeter is already present
    if ( !deepSpaceNetworkAdded_ )
    {
        // Altimeter has been created
        deepSpaceNetworkAdded_ = true;

        // Make sure that Earth is present in the simulation
        if ( bodyMap_.count( "Earth" ) == 0 )
        {
            throw std::runtime_error( "Error in creation of DSN system for body " + spacecraftName_ +
                                      ". Earth is not present in the simulated bodies, thus no tracking can be performed." );
        }

        // Generate random noise distribution
        generateDeepSpaceNetworkRandomNoiseDistribution( positionAccuracy, velocityAccuracy, lightTimeAccuracy );

        // Create function for computing corrupted spacecraft orientation
        deepSpaceNetworkFunction_ = boost::bind( &InstrumentsModel::getCurrentDeepSpaceNetworkTracking, this );
    }
    else
    {
        throw std::runtime_error( "Error in creation of DSN system for body " + spacecraftName_ +
                                  ". DSN tracking is already active." );
    }
}

//! Function to add a generic ranging system.
void InstrumentsModel::addGenericRangingSystem( const Eigen::Vector3d& positionBias,
                                                const Eigen::Vector3d& positionScaleFactor,
                                                const Eigen::Vector6d& positionMisalignment,
                                                const Eigen::Vector3d& positionAccuracy )
{
    // Check whether a generic ranging system is already present
    if ( !genericRangingSystemAdded_ )
    {
        // Inertial measurement unit has been created
        genericRangingSystemAdded_ = true;

        // Generate random noise distribution
        generateGenericRangingSystemRandomNoiseDistribution( positionAccuracy );

        // Create function for computing corrupted translational accelerations
        genericRangingSystemFunction_ = boost::bind(
                    &InstrumentsModel::getCurrentPosition, this, positionBias,
                    computeScaleMisalignmentMatrix( positionScaleFactor, positionMisalignment ) );
    }
    else
    {
        throw std::runtime_error( "Error in creation of inertial measurement unit for body " + spacecraftName_ +
                                  ". An IMU is already present." );
    }
}

//! Function to generate the noise distributions for the inertial measurement unit.
void InstrumentsModel::generateInertialMeasurementUnitRandomNoiseDistribution( const Eigen::Vector3d& accelerometerAccuracy,
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
                            normal_boost_distribution, { 0.0, accelerometerAccuracy[ i ] }, i ) );
        }
        else
        {
            accelerometerNoiseDistribution_.push_back( nullptr );
        }
    }

    // Create gyroscope noise distribution
    for ( unsigned int i = 0; i < 3; i++ )
    {
        if ( gyroscopeAccuracy[ i ] != 0.0 )
        {
            gyroscopeNoiseDistribution_.push_back(
                        createBoostContinuousRandomVariableGenerator(
                            normal_boost_distribution, { 0.0, gyroscopeAccuracy[ i ] }, 3 + i ) );
        }
        else
        {
            gyroscopeNoiseDistribution_.push_back( nullptr );
        }
    }
}

//! Function to generate the noise distributions for the star trackers.
void InstrumentsModel::generateStarTrackerRandomNoiseDistribution( const Eigen::Vector3d& starTrackerAccuracy )
{
    using namespace tudat::statistics;

    // Create star tracker noise distribution
    for ( unsigned int i = 0; i < 3; i++ )
    {
        if ( starTrackerAccuracy[ i ] != 0.0 )
        {
            starTrackerNoiseDistribution_.push_back(
                        createBoostContinuousRandomVariableGenerator(
                            normal_boost_distribution, { 0.0, starTrackerAccuracy[ i ] }, 6 + i ) );
        }
        else
        {
            starTrackerNoiseDistribution_.push_back( nullptr );
        }
    }
}

//! Function to generate the noise distributions for the altimeter.
void InstrumentsModel::generateAltimeterRandomNoiseDistribution( const double altimeterAccuracy )
{
    using namespace tudat::statistics;

    // Create altimeter noise distribution
    if ( altimeterAccuracy != 0.0 )
    {
        altimeterNoiseDistribution_ = createBoostContinuousRandomVariableGenerator(
                    normal_boost_distribution, { 0.0, altimeterAccuracy }, 9 );
    }
    else
    {
        altimeterNoiseDistribution_ = nullptr;
    }
}

//! Function to generate the noise distributions for the Deep Space Network system.
void InstrumentsModel::generateDeepSpaceNetworkRandomNoiseDistribution( const double positionAccuracy,
                                                                        const double velocityAccuracy,
                                                                        const double lightTimeAccuracy )
{
    using namespace tudat::statistics;

    // Combine accuracies
    Eigen::Vector2d combinedAccuracies;
    combinedAccuracies[ 0 ] = positionAccuracy;
    combinedAccuracies[ 1 ] = velocityAccuracy;
    combinedAccuracies[ 2 ] = lightTimeAccuracy;

    // Create Deep Space Network noise distribution
    for ( unsigned int i = 0; i < 3; i++ )
    {
        if ( combinedAccuracies[ i ] != 0.0 )
        {
            deepSpaceNetworkNoiseDistribution_.push_back(
                        createBoostContinuousRandomVariableGenerator(
                            normal_boost_distribution, { 0.0, combinedAccuracies[ i ] }, 10 + i ) );
        }
        else
        {
            deepSpaceNetworkNoiseDistribution_.push_back( nullptr );
        }
    }
}

//! Function to generate the noise distributions for the generic ranging system.
void InstrumentsModel::generateGenericRangingSystemRandomNoiseDistribution( const Eigen::Vector3d& positionAccuracy )
{
    using namespace tudat::statistics;

    // Create star tracker noise distribution
    for ( unsigned int i = 0; i < 3; i++ )
    {
        if ( positionAccuracy[ i ] != 0.0 )
        {
            positionNoiseDistribution_.push_back(
                        createBoostContinuousRandomVariableGenerator(
                            normal_boost_distribution, { 0.0, positionAccuracy[ i ] }, 13 + i ) );
        }
        else
        {
            positionNoiseDistribution_.push_back( nullptr );
        }
    }
}

} // namespace system_models

} // namespace tudat
