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

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"

namespace tudat
{

namespace system_models
{

//! Function to compute the scale-misalignment matrix from the scale and misalignment vectors.
Eigen::Matrix3d computeScaleMisalignmentMatrix( const Eigen::Vector3d& scaleFactorVector,
                                                const Eigen::Vector6d& misalignmentVector );

//! Class for guidance system of an aerobraking maneuver.
class NavigationInstrumentsModel
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param bodyMap
     *  \param selectedAccelerationPerBody
     *  \param centralBodies
     *  \param spacecraftName
     */
    NavigationInstrumentsModel( const simulation_setup::NamedBodyMap& bodyMap,
                                const simulation_setup::SelectedAccelerationMap& selectedAccelerationPerBody,
                                const std::map< std::string, std::string >& centralBodies,
                                const std::string& spacecraftName ) :
        bodyMap_( bodyMap ), selectedAccelerationPerBody_( selectedAccelerationPerBody ),
        centralBodies_( centralBodies ), spacecraftName_( spacecraftName ), currentTime_( 0.0 )
    {
        // Set instrument presence to false
        inertialMeasurementUnitAdded_ = false;
        starTrackerAdded_ = false;
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

    //! Function to update the onboard instruments to the current time.
    /*!
     *  Function to update the onboard instruments to the current time.
     *  \param currentTime Time at which the measurement needs to be taken.
     */
    void updateInstruments( const double currentTime )
    {
        // If current time has not been already updated
        if ( currentTime_ != currentTime )
        {
            // Update current time
            currentTime_ = currentTime;

            // Update inertial measurement unit
            if ( inertialMeasurementUnitAdded_ )
            {
                // Translational accelerations
                currentTranslationalAcceleration_.setZero( );
                inertialMeasurementUnitTranslationalAccelerationFunction_( currentTranslationalAcceleration_ );

                // Rotational velocity
                currentRotationalVelocity_.setZero( );
                inertialMeasurementUnitRotationalVelocityFunction_( currentRotationalVelocity_ );

                // Save inertial measurement unit measurements
                currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_[ currentTime_ ] =
                        ( Eigen::Vector6d << currentTranslationalAcceleration_, currentRotationalVelocity_ ).finished( );
            }

            // Update star tracker
            if ( starTrackerAdded_ )
            {
                // Compute and output orientation
                currentQuaternionToBaseFrame_.setZero( );
                starTrackerOrientationFunction_( currentQuaternionToBaseFrame_ );
                currentOrbitHistoryOfStarTrackerMeasurements_[ currentTime_ ] = currentQuaternionToBaseFrame_;
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
        return ( Eigen::Vector6d << getCurrentAccelerometerMeasurement( ),
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
            throw std::runtime_error( "Error while retrieving accelerations from onboard instrument system. No "
                                      "inertial measurement unit is present." );
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
            throw std::runtime_error( "Error while retrieving accelerations from onboard instrument system. No "
                                      "inertial measurement unit is present." );
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
        // Check that an inertial measurement unit is present in the spacecraft
        if ( starTrackerAdded_ )
        {
            return currentQuaternionToBaseFrame_;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving rotational velocity from onboard instrument system. No "
                                      "star trackers are present." );
        }
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
    }

private:

    //! Function to retrieve current translational accelerations exerted on the spacecraft.
    void getCurrentTranslationalAcceleration( Eigen::Vector3d& currentTranslationalAcceleration,
                                              const Eigen::Vector3d& biasVector, const Eigen::Matrix3d& scaleMisalignmentMatrix )
    {
        // Iterate over all accelerations acting on body
        basic_astrodynamics::SingleBodyAccelerationMap accelerationsOnBody = accelerationModelMap_.at( spacecraftName_ );
        for ( accelerationIterator_ = accelerationsOnBody.begin( ); accelerationIterator_ != accelerationsOnBody.end( );
              accelerationIterator_++ )
        {
            // Loop over each acceleration and disregard the central gravitational accelerations,
            // since IMUs do not measure them
            for ( unsigned int i = 1; i < accelerationIterator_->second.size( ); i++ )
            {
                // Calculate acceleration and add to state derivative
                currentTranslationalAcceleration += accelerationIterator_->second[ i ]->getAcceleration( );
            }
        }

        // Add errors to acceleration value
        currentTranslationalAcceleration = scaleMisalignmentMatrix * currentTranslationalAcceleration;
        currentTranslationalAcceleration += biasVector + produceAccelerometerNoise( );
    }

    //! Function to retrieve current rotational velocity of the spacecraft.
    void getCurrentRotationalVelocity( Eigen::Vector3d& currentRotationalVelocity,
                                       const Eigen::Vector3d& biasVector, const Eigen::Matrix3d& scaleMisalignmentMatrix )
    {
        // Iterate over all accelerations acting on body
        currentRotationalVelocity = bodyMap_.at( spacecraftName_ )->getCurrentAngularVelocityVectorInGlobalFrame( );

        // Add errors to acceleration value
        currentRotationalVelocity = scaleMisalignmentMatrix * currentRotationalVelocity;
        currentRotationalVelocity += biasVector + produceGyroscopeNoise( );
    }

    //! Function to retrieve current inertial orientation of the spacecraft.
    void getCurrentAttitude( Eigen::Vector4d& currentQuaternionToBaseFrame )
    {
        // Iterate over all accelerations acting on body
        currentQuaternionToBaseFrame = bodyMap_.at( spacecraftName_ )->getCurrentRotationToGlobalFrame( );

        // Add errors to acceleration value
        currentQuaternionToBaseFrame += produceStarTrackerNoise( );
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
    void generateStarTrackerRandomNoiseDistribution( const Eigen::Vector4d& starTrackerAccuracy );

    //! Function to produce accelerometer noise.
    /*!
     *  Function to produce accelerometer noise.
     *  \return Vector where the noise for the accelerometer of the inertial measurement unit is stored.
     */
    Eigen::Vector3d produceAccelerometerNoise( )
    {
        // Declare noise vector
        Eigen::Vector3d accelerometerNoise;

        // Loop over dimensions and add noise
        for ( int i = 0; i < 3; i++ )
        {
            if ( accelerometerNoiseDistribution_.at( i ) != NULL )
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
        Eigen::Vector3d gyroscopeNoise;

        // Loop over dimensions and add noise
        for ( int i = 0; i < 3; i++ )
        {
            if ( gyroscopeNoiseDistribution_.at( i ) != NULL )
            {
                gyroscopeNoise.second[ i ] = gyroscopeNoiseDistribution_.at( i )->getRandomVariableValue( );
            }
        }

        // Give back noise
        //        gyroscopeNoiseHistory_.push_back( gyroscopeNoise );
        return gyroscopeNoise;
    }

    //! Function to produce star tracker noise.
    /*!
     *  Function to produce star tracker noise.
     *  \return Vector where the noise for the star tracker is stored.
     */
    Eigen::Vector4d produceStarTrackerNoise( )
    {
        // Declare noise vector
        Eigen::Vector4d starTrackerNoise;

        // Loop over dimensions and add noise
        for ( int i = 0; i < 4; i++ )
        {
            if ( starTrackerNoiseDistribution_.at( i ) != NULL )
            {
                starTrackerNoise[ i ] = starTrackerNoiseDistribution_.at( i )->getRandomVariableValue( );
            }
        }

        // Give back noise
        //        starTrackerNoiseHistory_.push_back( starTrackerNoise );
        return starTrackerNoise;
    }

    //! Body map of the simulation.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Translational accelerations acting on the spacecraft.
    simulation_setup::SelectedAccelerationMap selectedAccelerationPerBody_;

    //! Central bodies of the simulation.
    std::map< std::string, std::string > centralBodies_;

    //! String denoting the name of the spacecraft body.
    std::string spacecraftName_;

    //! Double denoting current time.
    double currentTime_;

    //! Boolean denoting whether an inertial measurement unit is present in the spacecraft.
    bool inertialMeasurementUnitAdded_;

    //! Boolean denoting whether a star tracker is present in the spacecraft.
    bool starTrackerAdded_;

    //! Vector where the noise generators for the translational accelerations are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > accelerometerNoiseDistribution_;

    //! Vector where the noise generators for the rotational velocity are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > gyroscopeNoiseDistribution_;

    //! Vector where the noise generators for the spacecraft attitude are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > starTrackerNoiseDistribution_;

    //! Pointer to accelerations exerted on the spacecraft.
    basic_astrodynamics::AccelerationMap accelerationModelMap_;

    //! Predefined iterator to save (de)allocation time.
    basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationIterator_;

    //! Vector denoting current translational accelerations as measured by the accelerometer.
    Eigen::Vector3d currentTranslationalAcceleration_;

    //! Vector denoting current rotational velocity as measured by the gyroscope.
    Eigen::Vector3d currentRotationalVelocity_;

    //! Vector denoting current quaternion to base frame as measured by the star tracker.
    Eigen::Vector4d currentQuaternionToBaseFrame_;

    //! Function to compute the translational acceleration measured by the inertial measurement unit.
    boost::function< void( Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Matrix3d& ) >
    inertialMeasurementUnitTranslationalAccelerationFunction_;

    //! Function to compute the rotational velocity measured by the inertial measurement unit.
    boost::function< void( Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Matrix3d& ) >
    inertialMeasurementUnitRotationalVelocityFunction_;

    //! Function to compute the rotational velocity measured by the inertial measurement unit.
    boost::function< void( Eigen::Vector4d&, const Eigen::Vector3d& ) >
    starTrackerOrientationFunction_;

    //! Map of translational accelerations and rotational velocities measured by the inertial measurment unit
    //! during the current orbit.
    std::map< double, Eigen::Vector6d > currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_;

    //! Map of attitude measured by the star tracker during the current orbit.
    std::map< double, Eigen::Vector4d > currentOrbitHistoryOfStarTrackerMeasurements_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_NAVIGATION_INSTRUMENTS_MODEL_H
