/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NAVIGATION_INSTRUMENTS_MODEL_H
#define TUDAT_NAVIGATION_INSTRUMENTS_MODEL_H

#include <cmath>
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace system_models
{

//! Function to compute the scale-misalignment matrix from the scale and misalignment vectors.
Eigen::Matrix3d computeScaleMisalignmentMatrix( const Eigen::Vector3d& scaleFactorVector,
                                                const Eigen::Vector6d& misalignmentVector );

//! Function to add a small angle uncertainty to a quaternion vector.
Eigen::Vector4d sumQuaternionUncertainty( const Eigen::Vector4d& quaternionVector,
                                          const Eigen::Vector3d& equivalentVectorUncertainty );

//! Class for instrument models of a planet-bound spacecraft.
class NavigationInstrumentsModel
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param bodyMap Body map used in the main simulation.
     *  \param accelerationModelMap Acceleration model map used in the main simulation.
     *  \param spacecraftName Name of the spacecraft whose instrument model needs to be created.
     *  \param planetName Name of the planet the spacecraft is orbiting.
     */
    NavigationInstrumentsModel( const simulation_setup::NamedBodyMap& bodyMap,
                                const basic_astrodynamics::AccelerationMap& accelerationModelMap,
                                const std::string& spacecraftName,
                                const std::string& planetName ) :
        bodyMap_( bodyMap ), accelerationModelMap_( accelerationModelMap ),
        spacecraftName_( spacecraftName ), planetName_( planetName ), currentTime_( 0.0 )
    {
        // Set instrument presence to false
        inertialMeasurementUnitAdded_ = false;
        starTrackerAdded_ = false;
        altimeterAdded_ = false;

        // Get index of central body acceleration (which is not measured by the IMUs)
        for ( accelerationMapIterator_ = accelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != accelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                if ( ( basic_astrodynamics::getAccelerationModelType( accelerationMapIterator_->second[ i ] ) ==
                       basic_astrodynamics::spherical_harmonic_gravity ) &&
                     ( accelerationMapIterator_->first == planetName_ ) )
                {
                    sphericalHarmonicsGravityIndex_ = i;
                }
                else if ( ( basic_astrodynamics::getAccelerationModelType( accelerationMapIterator_->second[ i ] ) ==
                            basic_astrodynamics::third_body_central_gravity ) &&
                          ( accelerationMapIterator_->first == "Sun" ) )
                {
                    thirdBodyGravityIndex_ = i;
                }
            }
        }
    }

    //! Destructor.
    ~NavigationInstrumentsModel( ) { }

    //! Function to add an inertial measurement unit to the spacecraft set of instruments.
    /*!
     *  Function to add an inertial measurement unit to the spacecraft set of instruments.
     *  \param accelerometerBias Bias of the accelerometer.
     *  \param accelerometerScaleFactor Scale factor of the accelerometer.
     *  \param accelerometerMisalignment Misalignment of the accelerometer.
     *  \param accelerometerAccuracy Accuracy of accelerometer along each axis (3 sigma).
     *  \param gyroscopeBias Bias of the gyroscope.
     *  \param gyroscopeScaleFactor Scale factor of the gyroscope.
     *  \param gyroscopeMisalignment Misalignment of the gyroscope.
     *  \param gyroscopeAccuracy Accuracy of gyroscope along each axis (3 sigma).
     */
    void addInertialMeasurementUnit( const Eigen::Vector3d& accelerometerBias,
                                     const Eigen::Vector3d& accelerometerScaleFactor,
                                     const Eigen::Vector6d& accelerometerMisalignment,
                                     const Eigen::Vector3d& accelerometerAccuracy,
                                     const Eigen::Vector3d& gyroscopeBias,
                                     const Eigen::Vector3d& gyroscopeScaleFactor,
                                     const Eigen::Vector6d& gyroscopeMisalignment,
                                     const Eigen::Vector3d& gyroscopeAccuracy );

    //! Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
    /*!
     *  Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
     *  \param numberOfStarTrackers Number of star trackers to add to the spacecraft.
     *  \param starTrackerAccuracy Accuracy of star tracker along each axis (3 sigma).
     */
    void addStarTracker( const unsigned int numberOfStarTrackers,
                         const Eigen::Vector3d& starTrackerAccuracy );

    //! Function to add an altimeter to the spacecraft set of instruments.
    /*!
     *  Function to add an altimeter to the spacecraft set of instruments.
     *  \param fixedBodyFramePointingDirection
     *  \param altitudeRange
     *  \param altimeterAccuracyAsAFunctionOfAltitude
     */
    void addAltimeter( const Eigen::Vector3d& fixedBodyFramePointingDirection,
                       const std::pair< double, double >& altitudeRange,
                       const boost::function< double( double ) >& altimeterAccuracyAsAFunctionOfAltitude );

    //! Function to update the onboard instruments to the current time.
    /*!
     *  Function to update the onboard instruments to the current time. Note that this function needs to be called before retrieving
     *  the measurements at the current time, since the functions to get the measurements do not update the instruments.
     *  \param currentTime Time at which the measurement needs to be taken.
     */
    void updateInstruments( const double currentTime )
    {
        // If instruments have not been already updated for the current time
        if ( currentTime_ != currentTime )
        {
            // Update current time
            currentTime_ = currentTime;

            // Update inertial measurement unit
            if ( inertialMeasurementUnitAdded_ )
            {
                // Translational accelerations and rotational velocity
                inertialMeasurementUnitTranslationalAccelerationFunction_( );
                inertialMeasurementUnitRotationalVelocityFunction_( );

                // Save inertial measurement unit measurements
                currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_[ currentTime_ ] =
                        ( Eigen::VectorXd( 6 ) << currentTranslationalAcceleration_, currentRotationalVelocity_ ).finished( );
            }

            // Update star tracker
            if ( starTrackerAdded_ )
            {
                // Update and save attitude measurement
                starTrackerOrientationFunction_( );
                currentOrbitHistoryOfStarTrackerMeasurements_[ currentTime_ ] = currentQuaternionToBaseFrame_;
            }

            // Update altimeter
            if ( altimeterAdded_ )
            {
                // Update and save altimeter measurement
                altimeterFunction_( );
                currentOrbitHistoryOfAltimeterMeasurements_[ currentTime_ ] = currentAltitude_;
            }
        }
    }

    //! Function to retrieve current inertial measurement unit measurement.
    /*!
     *  Function to retrieve current inertial measurement unit measurement. This function simply calles
     *  the getCurrentAccelerometerMeasurement and getCurrentGyroscopeMeasurement functions to retireve
     *  the two IMU measurements separately and then merges them into one.
     *  \return Six-dimensional vector representing the translational accelerations and rotational
     *      velocities measured by the modeled inertial measurement unit.
     */
    Eigen::Vector6d getCurrentInertialMeasurementUnitMeasurement( )
    {
        // Merge translational and rotational accelerations and give result
        return ( Eigen::VectorXd( 6 ) << getCurrentAccelerometerMeasurement( ),
                 getCurrentGyroscopeMeasurement( ) ).finished( );
    }

    //! Function to retrieve only the current accelerometer measurement.
    /*!
     *  Function to retrieve only the current accelerometer measurement. The translational acceleration is
     *  retireved from the acceleration model and before being returned, it is corrupted with the accelerometer
     *  errors provided as an input to the addInertialMeasurementUnit function.
     *  \return Three-dimensional vector representing the translational accelerations measured by the
     *      modeled inertial measurement unit.
     */
    Eigen::Vector3d getCurrentAccelerometerMeasurement( )
    {
        // Check that an inertial measurement unit is present in the spacecraft
        if ( inertialMeasurementUnitAdded_ )
        {
            return currentTranslationalAcceleration_;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving translational accelerations from onboard instrument "
                                      "system. No inertial measurement unit is present." );
        }
    }

    //! Function to retrieve only the current gyroscope measurement.
    /*!
     *  Function to retrieve only the current gyroscope measurement. The rotational velocity is retireved
     *  from the body model and before being returned, it is corrupted with the gyroscope
     *  errors provided as an input to the addInertialMeasurementUnit function.
     *  \return Three-dimensional vector representing the rotational velocities measured by the
     *      modeled inertial measurement unit.
     */
    Eigen::Vector3d getCurrentGyroscopeMeasurement( )
    {
        // Check that an inertial measurement unit is present in the spacecraft
        if ( inertialMeasurementUnitAdded_ )
        {
            return currentRotationalVelocity_;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving rotational velocities from onboard instrument system. "
                                      "No inertial measurement unit is present." );
        }
    }

    //! Function to retrieve current star tracker measurement.
    /*!
     *  Function to retrieve current star tracker measurement. The attitude is retireved from the body
     *  model and before being returned, it is corrupted with the star tracker errors provided as an
     *  input to the addStarTracker function.
     *  \return Four-dimensional vector representing the current quaternion to base frame measured
     *      by the modeled star tracker, i.e., the rotation from body-fixed to inertial frame.
     */
    Eigen::Vector4d getCurrentStarTrackerMeasurement( )
    {
        // Check that a star tracker is present in the spacecraft
        if ( starTrackerAdded_ )
        {
            return currentQuaternionToBaseFrame_;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving attitude from onboard instrument system. No "
                                      "star trackers are present." );
        }
    }

    //! Function to retrieve current altimeter measurement.
    /*!
     *  Function to retrieve current altimeter measurement.
     *  \return
     */
    double getCurrentAltimeterMeasurement( )
    {
        // Check that an altimeter is present in the spacecraft
        if ( altimeterAdded_ )
        {
            return currentAltitude_;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving altitude from onboard instrument system. No "
                                      "altimeter is present." );
        }
    }

    //! Get standard deviation of the norm of accelerations measured by the accelerometer.
    double getStandardDeviationOfNormOfAccelerometerMeasurements( )
    {
        // Retrieve norm of accelerometer measurements
        std::vector< double > vectorOfNormOfAccelerometerMeasurements;
        for ( inertialMeasurementUnitMeasurementIterator_ = currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.begin( );
              inertialMeasurementUnitMeasurementIterator_ != currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.end( );
              inertialMeasurementUnitMeasurementIterator_++ )
        {
            vectorOfNormOfAccelerometerMeasurements.push_back(
                        inertialMeasurementUnitMeasurementIterator_->second.segment( 0, 3 ).norm( ) );
        }

        // Compute and output standard deviation
        return statistics::computeSampleVariance( vectorOfNormOfAccelerometerMeasurements );
    }

    //! Get history of inertial measurement unit measurements for the current orbit.
    std::map< double, Eigen::Vector6d > getCurrentOrbitHistoryOfInertialMeasurmentUnitMeasurements( )
    {
        return currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_;
    }

    //! Get history of star tracker measurements for the current orbit.
    std::map< double, Eigen::Vector4d > getCurrentOrbitHistoryOfStarTrackerMeasurements( )
    {
        return currentOrbitHistoryOfStarTrackerMeasurements_;
    }

    //! Get history of altimeter measurements for the current orbit.
    std::map< double, double > getCurrentOrbitHistoryOfAltimeterMeasurements( )
    {
        return currentOrbitHistoryOfAltimeterMeasurements_;
    }

    //! Clear histories of inertial measurmenet and star tracker measurements for current orbit.
    void clearCurrentOrbitMeasurementHistories( )
    {
        if ( inertialMeasurementUnitAdded_ )
        {
            currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.clear( );
        }
        if ( starTrackerAdded_ )
        {
            currentOrbitHistoryOfStarTrackerMeasurements_.clear( );
        }
        if ( altimeterAdded_ )
        {
            currentOrbitHistoryOfAltimeterMeasurements_.clear( );
        }
    }

private:

    //! Function to retrieve current translational accelerations exerted on the spacecraft.
    void getCurrentTranslationalAcceleration( const Eigen::Vector3d& biasVector, const Eigen::Matrix3d& scaleMisalignmentMatrix )
    {
        // Clear translational accelerations vector
        currentTranslationalAcceleration_.setZero( );

        // Iterate over all accelerations acting on body
        for ( accelerationMapIterator_ = accelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != accelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                // Disregard the central gravitational accelerations, since IMUs do not measure them
                if ( !( ( i == sphericalHarmonicsGravityIndex_ ) && ( accelerationMapIterator_->first == planetName_ ) ) &&
                     !( ( i == thirdBodyGravityIndex_ ) && ( accelerationMapIterator_->first == "Sun" ) ) )
                {
                    // Calculate acceleration and add to state derivative
                    currentTranslationalAcceleration_ += accelerationMapIterator_->second[ i ]->getAcceleration( );
                }
            }
        }

        // Add errors to acceleration value
        currentTranslationalAcceleration_ = scaleMisalignmentMatrix * currentTranslationalAcceleration_;
        currentTranslationalAcceleration_ += biasVector + produceAccelerometerNoise( );
    }

    //! Function to retrieve current rotational velocity of the spacecraft.
    void getCurrentRotationalVelocity( const Eigen::Vector3d& biasVector, const Eigen::Matrix3d& scaleMisalignmentMatrix )
    {
        // Iterate over all accelerations acting on body
        currentRotationalVelocity_ = bodyMap_.at( spacecraftName_ )->getCurrentAngularVelocityVectorInLocalFrame( );

        // Add errors to acceleration value
        currentRotationalVelocity_ = scaleMisalignmentMatrix * currentRotationalVelocity_;
        currentRotationalVelocity_ += biasVector + produceGyroscopeNoise( );
    }

    //! Function to retrieve current inertial orientation of the spacecraft.
    void getCurrentAttitude( )
    {
        // Get current attitude of spacecraft
        currentQuaternionToBaseFrame_ = linear_algebra::convertQuaternionToVectorFormat(
                    bodyMap_.at( spacecraftName_ )->getCurrentRotationToGlobalFrame( ) );

        // Add errors to attitude value
        sumQuaternionUncertainty( currentQuaternionToBaseFrame_, produceStarTrackerNoise( ) );
    }

    //! Function to retrieve current altitude of the spacecraft.
    void getCurrentAltitude( const Eigen::Vector3d& fixedBodyFramePointingDirection,
                             const std::pair< double, double >& altitudeRange,
                             const boost::function< double( double ) >& altimeterAccuracyAsAFunctionOfAltitude )
    {
        // Get current altitude of spacecraft
        Eigen::Vector3d currentRadialVector = bodyMap_.at( spacecraftName_ )->getState( ).segment( 0, 3 ) -
                bodyMap_.at( planetName_ )->getState( ).segment( 0, 3 );
        double currentRadialDistance = currentRadialVector.norm( );
        double planetaryRadius = bodyMap_.at( planetName_ )->getShapeModel( )->getAverageRadius( );
        double currentPointingAngle = 0.0;

        // Check that pointing angle is not too large, i.e., if Mars is still in sight
        double maximumPointingAngle = std::asin( planetaryRadius / currentRadialDistance );
        if ( currentPointingAngle < maximumPointingAngle )
        {
            // Get actual altitude measurement (i.e., accounting for the pointing angle)
            if ( std::fabs( currentPointingAngle ) < std::numeric_limits< double >::epsilon( ) )
            {
                currentAltitude_ = currentRadialDistance - planetaryRadius;
            }
            else
            {
                double sineOfCurrentPointingAngle = std::sin( currentPointingAngle );
                currentAltitude_ = planetaryRadius / sineOfCurrentPointingAngle *
                        std::sin( std::asin( currentRadialDistance / planetaryRadius * sineOfCurrentPointingAngle ) -
                                  currentPointingAngle );
            }
        }
        else
        {
            // Altimeter does not have Mars in sight
            currentAltitude_ = TUDAT_NAN;
        }

        // Check that altitude is within altimeter limits
        if ( ( currentAltitude_ < altitudeRange.first ) || ( currentAltitude_ > altitudeRange.second ) )
        {
            // Altitude cannot be measured
            currentAltitude_ = TUDAT_NAN; // replace with function that gives large errors?
        }

        // Add errors to altitude value
        currentAltitude_ += produceAltimeterNoise( ) * altimeterAccuracyAsAFunctionOfAltitude( currentAltitude_ );
    }

    //! Function to generate the noise distributions for the inertial measurement unit.
    /*!
     *  Function to generate the noise distributions for both instruments of the inertial measurement unit,
     *  which uses a Gaussian distribution, with zero mean and standard deviation given by the accuracy of the
     *  accelerometer and gyroscopes.
     *  \param accelerometerAccuracy Accuracy of accelerometer along each axis (3 sigma).
     *  \param gyroscopeAccuracy Accuracy of gyroscope along each axis (3 sigma).
     */
    void generateInertialMeasurementUnitRandomNoiseDistribution( const Eigen::Vector3d& accelerometerAccuracy,
                                                                 const Eigen::Vector3d& gyroscopeAccuracy );

    //! Function to generate the noise distributions for the star trackers.
    /*!
     *  Function to generate the noise distributions for the star trackers, which uses a Gaussian distribution,
     *  with zero mean and standard deviation given by the accuracy of the star tracker.
     *  \param starTrackerAccuracy Accuracy of star tracker along each axis (3 sigma).
     */
    void generateStarTrackerRandomNoiseDistribution( const Eigen::Vector3d& starTrackerAccuracy );

    //! Function to generate the noise distributions for the altimeter.
    /*!
     *  Function to generate the noise distributions for the altimeter, which uses a Gaussian distribution,
     *  with zero mean and standard deviation given by the accuracy of the altimeter.
     *  \param altimeterAccuracy Accuracy of altimeter (3 sigma).
     */
    void generateAltimeterRandomNoiseDistribution( const double altimeterAccuracy );

    //! Function to produce accelerometer noise.
    /*!
     *  Function to produce accelerometer noise.
     *  \return Vector where the noise for the accelerometer of the inertial measurement unit is stored.
     */
    Eigen::Vector3d produceAccelerometerNoise( )
    {
        // Declare noise vector
        Eigen::Vector3d accelerometerNoise = Eigen::Vector3d::Zero( );

        // Loop over dimensions and add noise
        for ( int i = 0; i < 3; i++ )
        {
            if ( accelerometerNoiseDistribution_.at( i ) != nullptr )
            {
                accelerometerNoise[ i ] = accelerometerNoiseDistribution_.at( i )->getRandomVariableValue( );
            }
        }

        // Give back noise
        //        accelerometerNoiseHistory_.push_back( accelerometerNoise );
        return accelerometerNoise;
    }

    //! Function to produce gyroscope noise.
    /*!
     *  Function to produce gyroscope noise.
     *  \return Vector where the noise for the gyroscope of the inertial measurement unit is stored.
     */
    Eigen::Vector3d produceGyroscopeNoise( )
    {
        // Declare noise vector
        Eigen::Vector3d gyroscopeNoise = Eigen::Vector3d::Zero( );

        // Loop over dimensions and add noise
        for ( int i = 0; i < 3; i++ )
        {
            if ( gyroscopeNoiseDistribution_.at( i ) != nullptr )
            {
                gyroscopeNoise[ i ] = gyroscopeNoiseDistribution_.at( i )->getRandomVariableValue( );
            }
        }

        // Give back noise
        //        gyroscopeNoiseHistory_.push_back( gyroscopeNoise );
        return gyroscopeNoise;
    }

    //! Function to produce star tracker noise.
    /*!
     *  Function to produce star tracker noise. Note that this noise is output in terms of angle-axis accuracy.
     *  Before being added (via quaternion multiplication) to the quaternion state, it has to be converted to
     *  quaternion. This is done via the sumQuaternionUncertainty function.
     *  \return Vector where the noise for the star tracker is stored.
     */
    Eigen::Vector3d produceStarTrackerNoise( )
    {
        // Declare noise vector
        Eigen::Vector3d starTrackerNoise = Eigen::Vector3d::Zero( );

        // Loop over dimensions and add noise
        for ( int i = 0; i < 3; i++ )
        {
            if ( starTrackerNoiseDistribution_.at( i ) != nullptr )
            {
                starTrackerNoise[ i ] = starTrackerNoiseDistribution_.at( i )->getRandomVariableValue( );
            }
        }

        // Give back noise
        //        starTrackerNoiseHistory_.push_back( starTrackerNoise );
        return starTrackerNoise;
    }

    //! Function to produce altimeter noise.
    /*!
     *  Function to produce altimeter noise.
     *  \return Double where the noise for the altimeter is stored.
     */
    double produceAltimeterNoise( )
    {
        // Declare noise value
        double altimeterNoise = 0;

        // Add noise
        if ( altimeterNoiseDistribution_ != nullptr )
        {
            altimeterNoise = altimeterNoiseDistribution_->getRandomVariableValue( );
        }

        // Give back noise
        //        altimeterNoiseHistory_.push_back( altimeterNoise );
        return altimeterNoise;
    }

    //! Body map of the simulation.
    const simulation_setup::NamedBodyMap bodyMap_;

    //! Pointer to accelerations exerted on the spacecraft.
    const basic_astrodynamics::AccelerationMap accelerationModelMap_;

    //! String denoting the name of the spacecraft body.
    const std::string spacecraftName_;

    //! String denoting the name of the planet being orbited body.
    const std::string planetName_;

    //! Double denoting current time.
    double currentTime_;

    //! Boolean denoting whether an inertial measurement unit is present in the spacecraft.
    bool inertialMeasurementUnitAdded_;

    //! Boolean denoting whether a star tracker is present in the spacecraft.
    bool starTrackerAdded_;

    //! Boolean denoting whether an altimeter is present in the spacecraft.
    bool altimeterAdded_;

    //! Integer denoting the index of the spherical harmonics gravity exerted by the planet being orbited.
    unsigned int sphericalHarmonicsGravityIndex_;

    //! Integer denoting the index of the third body gravity exerted by the Sun.
    unsigned int thirdBodyGravityIndex_;

    //! Vector where the noise generators for the translational accelerations measurement are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > accelerometerNoiseDistribution_;

    //! Vector where the noise generators for the rotational velocity measurement are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > gyroscopeNoiseDistribution_;

    //! Vector where the noise generators for the spacecraft attitude measurement are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > starTrackerNoiseDistribution_;

    //! Noise generator for the spacecraft altitude measurement.
    boost::shared_ptr< statistics::RandomVariableGenerator< double > > altimeterNoiseDistribution_;

    //! Vector denoting current translational accelerations as measured by the accelerometer.
    Eigen::Vector3d currentTranslationalAcceleration_;

    //! Vector denoting current rotational velocity as measured by the gyroscope.
    Eigen::Vector3d currentRotationalVelocity_;

    //! Vector denoting current quaternion to base frame as measured by the star tracker.
    Eigen::Vector4d currentQuaternionToBaseFrame_;

    //! Vector denoting current altitude as measured by the altimeter.
    double currentAltitude_;

    //! Function to compute the translational acceleration measured by the inertial measurement unit.
    boost::function< void( ) > inertialMeasurementUnitTranslationalAccelerationFunction_;

    //! Function to compute the rotational velocity measured by the inertial measurement unit.
    boost::function< void( ) > inertialMeasurementUnitRotationalVelocityFunction_;

    //! Function to compute the attitude measured by the star tracker.
    boost::function< void( ) > starTrackerOrientationFunction_;

    //! Function to compute the altitude measured by the altimeter.
    boost::function< void( ) > altimeterFunction_;

    //! Map of translational accelerations and rotational velocities measured by the inertial measurment unit
    //! during the current orbit.
    std::map< double, Eigen::Vector6d > currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_;

    //! Map of attitude measured by the star tracker during the current orbit.
    std::map< double, Eigen::Vector4d > currentOrbitHistoryOfStarTrackerMeasurements_;

    //! Map of attitude measured by the star tracker during the current orbit.
    std::map< double, double > currentOrbitHistoryOfAltimeterMeasurements_;

    //! Predefined iterator to save (de)allocation time.
    basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationMapIterator_;

    //! Predefined iterator to save (de)allocation time.
    std::map< double, Eigen::Vector6d >::const_iterator inertialMeasurementUnitMeasurementIterator_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_NAVIGATION_INSTRUMENTS_MODEL_H
