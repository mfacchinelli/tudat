/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NAVIGATION_INSTRUMENT_H
#define TUDAT_NAVIGATION_INSTRUMENT_H

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

//! Class for guidance system of an aerobraking maneuver.
class NavigationInstrumentModel
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
    NavigationInstrumentModel( const simulation_setup::NamedBodyMap& bodyMap,
                               const simulation_setup::SelectedAccelerationMap& selectedAccelerationPerBody,
                               const std::map< std::string, std::string >& centralBodies,
                               const std::string& spacecraftName ) :
        bodyMap_( bodyMap ), selectedAccelerationPerBody_( selectedAccelerationPerBody ),
        centralBodies_( centralBodies ), spacecraftName_( spacecraftName )
    {
        // Set instrument presence to false
        inertialMeasurementUnitAdded_ = false;
        starTrackerAdded_ = false;
    }

    //! Destructor.
    ~NavigationInstrumentModel( ) { }

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
                                     const Eigen::Vector3d& gyroscopeAccuracy )
    {
        // Check whether an inertial measurement unit is already present
        if ( !inertialMeasurementUnitAdded_ )
        {
            // Inertial measurement unit has been created
            inertialMeasurementUnitAdded_ = true;

            // Generate random noise distribution
            generateInertialMeasurementUnitRandomNoiseDistribution( accelerometerAccuracy, gyroscopeAccuracy );

            // Create acceleration model object
            accelerationModelMap_ = simulation_setup::createAccelerationModelsMap( bodyMap_, selectedAccelerationPerBody_,
                                                                                   centralBodies_ );

            // Create function for computing corrupted translational accelerations
            inertialMeasurementUnitTranslationalAccelerationFunction_ = boost::bind(
                        &NavigationInstrumentModel::getCurrentTranslationalAcceleration, this, _1, accelerometerBias,
                        computeScaleMisalignmentMatrix( accelerometerScaleFactor, accelerometerMisalignment ) );

            // Create function for computing corrupted rotational velocity
            inertialMeasurementUnitRotationalVelocityFunction_ = boost::bind(
                        &NavigationInstrumentModel::getCurrentRotationalVelocity, this, _1, gyroscopeBias,
                        computeScaleMisalignmentMatrix( gyroscopeScaleFactor, gyroscopeMisalignment ) );
        }
        else
        {
            throw std::runtime_error( "Error in creation of inertial measurement unit for body " + spacecraftName_ +
                                      ". An IMU is already present." );
        }
    }

    //! Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
    /*!
     *  Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
     *  \param numberOfStarTrackers Number of star trackers to add to the spacecraft.
     *  \param starTrackerAccuracy Accuracy of star tracker along each axis (3 sigma).
     */
    void addStarTracker( const unsigned int numberOfStarTrackers,
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
                            &NavigationInstrumentModel::getCurrentAttitude, this, _1 );
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

    //! Function to retrieve current inertial measurement unit measurement.
    Eigen::Vector6d getCurrentInertialMeasurementUnitMeasurement( )
    {
        // Check that an inertial measurement unit is present in the spacecraft
        if ( inertialMeasurementUnitAdded_ )
        {
            // Define output vector
            Eigen::Vector6d currentInertialMeasurementUnitMeasurement;

            // Translational accelerations
            Eigen::Vector3d currentTranslationalAcceleration = Eigen::Vector3d::Zero( );
            inertialMeasurementUnitTranslationalAccelerationFunction_( currentTranslationalAcceleration );

            // Rotational velocity
            Eigen::Vector3d currentRotationalVelocity;
            inertialMeasurementUnitRotationalVelocityFunction_( currentRotationalVelocity );

            // Merge translational and rotational accelerations
            currentInertialMeasurementUnitMeasurement << currentTranslationalAcceleration, currentRotationalVelocity;
            return currentInertialMeasurementUnitMeasurement;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving accelerations from onboard instrument system. No "
                                      "inertial measurement unit is present." );
        }
    }

    //! Function to retrieve current star tracker measurement.
    Eigen::Vector4d getCurrentStarTrackerMeasurement( )
    {
        // Check that an inertial measurement unit is present in the spacecraft
        if ( starTrackerAdded_ )
        {
            // Compute and output orientation
            Eigen::Vector4d currentQuaternionToBaseFrame;
            starTrackerOrientationFunction_( currentQuaternionToBaseFrame );
            return currentQuaternionToBaseFrame;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving rotational velocity from onboard instrument system. No "
                                      "star trackers are present." );
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
    /*!
     *  Function to generate the noise distributions for the star trackers, which uses a Gaussian distribution,
     *  with zero mean and standard deviation given by the accuracy of the star tracker.
     *  \param starTrackerAccuracy Accuracy of star tracker along each axis (3 sigma).
     */
    void generateInertialMeasurementUnitRandomNoiseDistribution( const Eigen::Vector4d& starTrackerAccuracy )
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

    //! Function to compute the translational acceleration measured by the inertial measurement unit.
    boost::function< void( Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Matrix3d& ) >
    inertialMeasurementUnitTranslationalAccelerationFunction_;

    //! Function to compute the rotational velocity measured by the inertial measurement unit.
    boost::function< void( Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Matrix3d& ) >
    inertialMeasurementUnitRotationalVelocityFunction_;

    //! Function to compute the rotational velocity measured by the inertial measurement unit.
    boost::function< void( Eigen::Vector4d&, const Eigen::Vector3d& ) >
    starTrackerOrientationFunction_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_NAVIGATION_INSTRUMENT_H
