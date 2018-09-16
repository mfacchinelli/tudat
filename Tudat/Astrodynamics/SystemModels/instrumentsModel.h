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

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"
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
class InstrumentsModel
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
    InstrumentsModel( const simulation_setup::NamedBodyMap& bodyMap,
                      const basic_astrodynamics::AccelerationMap& accelerationModelMap,
                      const std::string& spacecraftName,
                      const std::string& planetName ) :
        bodyMap_( bodyMap ), accelerationModelMap_( accelerationModelMap ),
        spacecraftName_( spacecraftName ), planetName_( planetName ), currentTime_( TUDAT_NAN )
    {
        // Set instrument presence to false
        inertialMeasurementUnitAdded_ = false;
        starTrackerAdded_ = false;
        altimeterAdded_ = false;
        deepSpaceNetworkAdded_ = false;

        // Get index of central body acceleration (which is not measured by the IMUs)
        sphericalHarmonicsGravityIndex_ = static_cast< unsigned int >( TUDAT_NAN );
        thirdBodyGravityIndex_ = static_cast< unsigned int >( TUDAT_NAN );
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
    ~InstrumentsModel( ) { }

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
     *  \param fixedBodyFramePointingDirection Direction in the body-fixed frame in which the altimeter is pointing.
     *  \param altitudeRange Range of altitudes between which the altimeter provides measurements. It is defined as a pair, such that
     *      the first entry gives the lower altitude bound and the second one the upper bound.
     *  \param altimeterAccuracyAsAFunctionOfAltitude Function returning the accuracy of the altimeter as a function of altitude.
     *      This can be a constant value if defined with boost::lambda::constant.
     */
    void addAltimeter( const Eigen::Vector3d& fixedBodyFramePointingDirection,
                       const std::pair< double, double >& altitudeRange,
                       const boost::function< double( const double ) >& altimeterAccuracyAsAFunctionOfAltitude );

    //! Function to add a Deep Space Network measurement system.
    /*!
     *  Function to add a Deep Space Network measurement system, which acts by retrieving the position of the spacecraft
     *  at a previous time, given by subtracting the light time error from the current time.
     *  \param positionAccuracy Accuracy in determining spacecraft position (3 sigma).
     *  \param velocityAccuracy Accuracy in determining spacecraft velocity (3 sigma).
     *  \param lightTimeAccuracy Accuracy in knowledge of position of Earth and Mars, expressed in seconds (could be compared to
     *      the clock error of GNSS).
     */
    void addDeepSpaceNetwork( const double positionAccuracy,
                              const double velocityAccuracy,
                              const double lightTimeAccuracy );

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

            // Update Deep Space Network
            if ( deepSpaceNetworkAdded_ )
            {
                // Update and save Deep Space Network measurement
                deepSpaceNetworkFunction_( );
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

    //! Function to retrieve only the accelerometer measurement smoothed over time.
    /*!
     *  Function to retrieve only the accelerometer measurement smoothed over time. The translational acceleration is
     *  retireved from the acceleration model and before being returned, it is corrupted with the accelerometer
     *  errors provided as an input to the addInertialMeasurementUnit function.
     *  \return Three-dimensional vector representing the smoothed translational accelerations measured by the
     *      modeled inertial measurement unit.
     */
    Eigen::Vector3d getSmoothedAccelerometerMeasurement( )
    {
        // Check that an inertial measurement unit is present in the spacecraft
        if ( inertialMeasurementUnitAdded_ )
        {
            // Compute smoothed accelerometer measurement
            unsigned int limitingValue = ( currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.size( ) < 500 ) ?
                        currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.size( ) : 500;
            Eigen::Vector3d smoothedAccelerometerMeasurement = Eigen::Vector3d::Zero( );
            for ( inertialMeasurementUnitMeasurementIterator_ =
                  std::prev( currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.end( ), limitingValue );
                  inertialMeasurementUnitMeasurementIterator_ != currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.end( );
                  inertialMeasurementUnitMeasurementIterator_++ )
            {
                smoothedAccelerometerMeasurement += inertialMeasurementUnitMeasurementIterator_->second.segment( 0, 3 );
            }
            smoothedAccelerometerMeasurement /= limitingValue;
            return smoothedAccelerometerMeasurement;
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

    //! Function to retrieve only the gyroscope measurement smoothed over time.
    /*!
     *  Function to retrieve only the gyroscope measurement smoothed over time. The rotational velocity is retireved
     *  from the body model and before being returned, it is corrupted with the gyroscope
     *  errors provided as an input to the addInertialMeasurementUnit function.
     *  \return Three-dimensional vector representing the smoothed rotational velocities measured by the
     *      modeled inertial measurement unit.
     */
    Eigen::Vector3d getSmoothedGyroscopeMeasurement( )
    {
        // Check that an inertial measurement unit is present in the spacecraft
        if ( inertialMeasurementUnitAdded_ )
        {
            // Compute smoothed gyroscope measurement
            unsigned int limitingValue = ( currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.size( ) < 500 ) ?
                        currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.size( ) : 500;
            Eigen::Vector3d smoothedGyroscopeMeasurement = Eigen::Vector3d::Zero( );
            for ( inertialMeasurementUnitMeasurementIterator_ =
                  std::prev( currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.end( ), limitingValue );
                  inertialMeasurementUnitMeasurementIterator_ != currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.end( );
                  inertialMeasurementUnitMeasurementIterator_++ )
            {
                smoothedGyroscopeMeasurement += inertialMeasurementUnitMeasurementIterator_->second.segment( 0, 3 );
            }
            smoothedGyroscopeMeasurement /= limitingValue;
            return smoothedGyroscopeMeasurement;
        }
        else
        {
            throw std::runtime_error( "Error while retrieving translational accelerations from onboard instrument "
                                      "system. No inertial measurement unit is present." );
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
     *  \return Double denoting the current altimeter measurement, affected by noise. Note that if the current altitude is outside
     *      the range defined by altitudeRange, the result will be nonsensical.
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

    //! Function to retrieve current Deep Space Network measurement.
    /*!
     *  Function to retrieve current Deep Space Network measurement, as measured by the antennas and post-processed by the
     *  ground segment on Earth. Note that the measurement represents the state of the spacecraft at the time of last contact with
     *  Earth (i.e., 2 times the light-time distance to Earth).
     *  \return Pair of double and vector, denoting the light-time delay at the time of measurement (2 times the light-time distance
     *      to the spacecraft), and the Cartesian state computed by merging the tracking result with a post-processing software for
     *      orbit determination (done on the ground).
     */
    std::pair< double, Eigen::Vector6d > getCurrentDeepSpaceNetworkMeasurement( )
    {
        // Check that a Deep Space Network tracking is active
        if ( deepSpaceNetworkAdded_ )
        {
            // Retrieve light-time and states over time
            std::vector< double > vectorOfTrackingTimes = utilities::createVectorFromMapKeys( historyOfDeepSpaceNetworkMeasurements_ );
            std::vector< std::pair< double, Eigen::Vector6d > > vectorOfTrackedElements =
                    utilities::createVectorFromMapValues( historyOfDeepSpaceNetworkMeasurements_ );
            std::vector< double > vectorOfLightTimeDistances;
            std::vector< Eigen::Vector6d > vectorOfTrackedStates;
            for ( unsigned int i = 0; i < vectorOfTrackedElements.size( ); i++ )
            {
                vectorOfLightTimeDistances.push_back( vectorOfTrackedElements.at( i ).first );
                vectorOfTrackedStates.push_back( vectorOfTrackedElements.at( i ).second );
            }

            // Create interpolator object based on stored Deep Space Network tracking
            double currentDeepSpaceNetworkTrackingTime;
            Eigen::Vector6d currentDeepSpaceNetworkTrackedState;
            try
            {
                currentDeepSpaceNetworkTrackingTime = interpolators::CubicSplineInterpolator< double, double >(
                            vectorOfTrackingTimes, vectorOfLightTimeDistances, interpolators::huntingAlgorithm,
                            interpolators::throw_exception_at_boundary ).interpolate( currentTime_ );
                currentDeepSpaceNetworkTrackedState = interpolators::CubicSplineInterpolator< double, Eigen::Vector6d >(
                            vectorOfTrackingTimes, vectorOfTrackedStates, interpolators::huntingAlgorithm,
                            interpolators::throw_exception_at_boundary ).interpolate( currentTime_ );
            }
            catch ( const std::exception& caughtException )
            {
                std::cerr << caughtException.what( ) << std::endl;
                throw std::runtime_error( "Error while retrieving position and velocity from DSN measurement. No measurement is "
                                          "available at the time requested. See error above for details on the boundary values." );
            }

            // Give output as pair
            return std::make_pair( currentDeepSpaceNetworkTrackingTime, currentDeepSpaceNetworkTrackedState );
        }
        else
        {
            throw std::runtime_error( "Error while retrieving position and velocity from DSN measurement. "
                                      "DSN tracking is not active." );
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

    //! Get history of altimeter measurements for the current orbit.
    std::map< double, double > getCurrentOrbitHistoryOfAltimeterMeasurements( )
    {
        return currentOrbitHistoryOfAltimeterMeasurements_;
    }

    //! Function to randomize the perturbation coefficients.
    /*!
     *  Function to randomize the perturbation coefficients, by accessing the atmosphere of the body, transforming it
     *  to tabulated atmosphere and randomizing the coefficient vector.
     */
    void randomizeAtmospherePerturbations( )
    {
        // Access planet atmosphere and convert it to tabulated atmosphere
        boost::shared_ptr< aerodynamics::TabulatedAtmosphere > atmosphere =
                boost::dynamic_pointer_cast< aerodynamics::TabulatedAtmosphere >( bodyMap_.at( planetName_ )->getAtmosphereModel( ) );

        // Randomize perturbations layer
        if ( atmosphere != nullptr )
        {
            atmosphere->randomizeAtmospherePerturbations( );
        }
        else
        {
            throw std::runtime_error( "Error in instruments model. The atmosphere model is not a tabulated atmosphere." );
        }
    }

    //! Clear histories of inertial measurmenet and star tracker measurements for current orbit.
    void clearCurrentOrbitMeasurementHistories( )
    {
        // Erase history of inertial measurement unit measurements
        if ( inertialMeasurementUnitAdded_ )
        {
            currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.clear( );
        }

        // Erase history of star tracker measurements
        if ( starTrackerAdded_ )
        {
            currentOrbitHistoryOfStarTrackerMeasurements_.clear( );
        }

        // Erase history of altimeter measurements
        if ( altimeterAdded_ )
        {
            currentOrbitHistoryOfAltimeterMeasurements_.clear( );
        }

        // Erase history of Deep Space Network measurements
        if ( deepSpaceNetworkAdded_ )
        {
            // Loop over measurements taken
            for ( std::map< double, std::pair< double, Eigen::Vector6d > >::const_iterator
                  measurementIterator = historyOfDeepSpaceNetworkMeasurements_.begin( );
                  measurementIterator != historyOfDeepSpaceNetworkMeasurements_.end( ); measurementIterator++ )
            {
                if ( measurementIterator->first >= ( currentTime_ - 0.0417 * physical_constants::JULIAN_DAY ) )
                {
                    // Erase measurements made more than approximately 1 hour before current time
                    historyOfDeepSpaceNetworkMeasurements_.erase( historyOfDeepSpaceNetworkMeasurements_.begin( ),
                                                                  measurementIterator-- );
                    break;
                }
            }
        }
    }

    //! Function to revert to the previous time step.
    /*!
     *  Function to revert to the previous time step. This function is run if the current propagation needs to be stopped, since
     *  the current time will be run the next time the GNC system is called.
     *  \param currentTime Double denoting the current time, i.e., the instant that has to be discarded.
     */
    void revertToPreviousTimeStep( const double currentTime )
    {
        // Reset time
        currentTime_ = TUDAT_NAN;

        // Erase measurement of inertial measurement unit
        if ( inertialMeasurementUnitAdded_ )
        {
            currentTranslationalAcceleration_.setConstant( TUDAT_NAN );
            currentRotationalVelocity_.setConstant( TUDAT_NAN );
            if ( currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.count( currentTime ) != 0 )
            {
                currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_.erase( currentTime );
            }
        }

        // Erase measurement of star tracker
        if ( starTrackerAdded_ )
        {
            currentQuaternionToBaseFrame_.setConstant( TUDAT_NAN );
            if ( currentOrbitHistoryOfStarTrackerMeasurements_.count( currentTime ) != 0 )
            {
                currentOrbitHistoryOfStarTrackerMeasurements_.erase( currentTime );
            }
        }

        // Erase measurement of altimeter
        if ( altimeterAdded_ )
        {
            currentAltitude_ = TUDAT_NAN;
            if ( currentOrbitHistoryOfAltimeterMeasurements_.count( currentTime ) != 0 )
            {
                currentOrbitHistoryOfAltimeterMeasurements_.erase( currentTime );
            }
        }

        // Erase measurement of Deep Space Network
        if ( deepSpaceNetworkAdded_ )
        {
            if ( historyOfDeepSpaceNetworkMeasurements_.count( currentTime ) != 0 )
            {
                historyOfDeepSpaceNetworkMeasurements_.erase( currentTime );
            }
        }
    }

    Eigen::Vector6d getActualSpacecraftTranslationalState( )
    {
        return bodyMap_.at( spacecraftName_ )->getState( ) - bodyMap_.at( planetName_ )->getState( );
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
        currentTranslationalAcceleration_.noalias( ) += biasVector + produceAccelerometerNoise( );
    }

    //! Function to retrieve current rotational velocity of the spacecraft.
    void getCurrentRotationalVelocity( const Eigen::Vector3d& biasVector, const Eigen::Matrix3d& scaleMisalignmentMatrix )
    {
        // Iterate over all accelerations acting on body
        currentRotationalVelocity_ = bodyMap_.at( spacecraftName_ )->getCurrentAngularVelocityVectorInLocalFrame( );

        // Add errors to acceleration value
        currentRotationalVelocity_ = scaleMisalignmentMatrix * currentRotationalVelocity_;
        currentRotationalVelocity_.noalias( ) += biasVector + produceGyroscopeNoise( );
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
        Eigen::Vector3d currentRadialVector = bodyMap_.at( spacecraftName_ )->getPosition( ) - bodyMap_.at( planetName_ )->getPosition( );
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

    //! Function to retrieve current position and velocity of the spacecraft.
    void getCurrentDeepSpaceNetworkTracking( )
    {
        // Retrieve noise for the current time
        Eigen::Vector7d currentNoiseVector = produceDeepSpaceNetworkNoise( );

        // Get current position and velocity of the spacecraft
        Eigen::Vector6d currentState = bodyMap_.at( spacecraftName_ )->getState( ) - bodyMap_.at( planetName_ )->getState( );

        // Compute current light-time distance to Earth
        double currentPlanetaryRange = ( bodyMap_.at( "Earth" )->getPosition( ) -
                                         bodyMap_.at( spacecraftName_ )->getPosition( ) ).norm( );
        double currentLightTimeDistance = currentPlanetaryRange / physical_constants::SPEED_OF_LIGHT;

        // Add noise to the state and time
        currentState.noalias( ) += currentNoiseVector.segment( 0, 6 );
        currentLightTimeDistance += currentNoiseVector[ 6 ];

        // Store value to map of measurements
        // Note that this method of storing the state is not entirely correct. The light-time distance (LTD) of the instant when the
        // measurement is received will be different, and should be accounted for. Here, the current LTD is multiplied by two, but
        // is going to produce a slight error. Since Mars is never further than approximately 20 minutes, the error is negligible.
        historyOfDeepSpaceNetworkMeasurements_[ currentTime_ + 2.0 * currentLightTimeDistance ] =
                std::make_pair( currentLightTimeDistance, currentState );
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

    //! Function to generate the noise distributions for the Deep Space Network system.
    /*!
     *  Function to generate the noise distributions for the Deep Space Network system, which uses a Gaussian distribution,
     *  with zero mean and standard deviation given by the input accuracies.
     *  \param positionAccuracy Accuracy in determining spacecraft position (3 sigma).
     *  \param velocityAccuracy Accuracy in determining spacecraft velocity (3 sigma).
     *  \param lightTimeAccuracy Accuracy in knowledge of position of Earth and Mars, expressed in seconds (could be compared to
     *      the clock error of GNSS).
     */
    void generateDeepSpaceNetworkRandomNoiseDistribution( const double positionAccuracy,
                                                          const double velocityAccuracy,
                                                          const double lightTimeAccuracy );

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
        for ( unsigned int i = 0; i < 3; i++ )
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
        for ( unsigned int i = 0; i < 3; i++ )
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
        for ( unsigned int i = 0; i < 3; i++ )
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

    //! Function to produce Deep Space Network noise.
    /*!
     *  Function to produce Deep Space Network noise, where the position and velocity accuracies are used three times each, to
     *  produce a vector of dimension 6, to match the Cartesian coordinate vector.
     *  \return Vector where the noise for the Deep Space Network is stored.
     */
    Eigen::Vector7d produceDeepSpaceNetworkNoise( )
    {
        // Declare noise value
        Eigen::Vector7d deepSpaceNetworkNoise = Eigen::Vector7d::Zero( );

        // Loop over dimensions and add noise
        for ( unsigned int i = 0; i < 3; i++ )
        {
            unsigned int limitingValue = ( i != 2 ) ? 3 : 1;
            for ( unsigned int j = 0; j < limitingValue; j++ )
            {
                if ( deepSpaceNetworkNoiseDistribution_.at( i ) != nullptr )
                {
                    deepSpaceNetworkNoise[ 3 * i + j ] = deepSpaceNetworkNoiseDistribution_.at( i )->getRandomVariableValue( );
                }
            }
        }

        // Give back noise
        //        deepSpaceNetworkNoiseHistory_.push_back( deepSpaceNetworkNoise );
        return deepSpaceNetworkNoise;
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

    //! Boolean denoting whether Deep Space Network tracking is active.
    bool deepSpaceNetworkAdded_;

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

    //! Vector where the noise generators for the Deep Space Network are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > deepSpaceNetworkNoiseDistribution_;

    //! Vector denoting current translational accelerations as measured by the accelerometer.
    Eigen::Vector3d currentTranslationalAcceleration_;

    //! Vector denoting current rotational velocity as measured by the gyroscope.
    Eigen::Vector3d currentRotationalVelocity_;

    //! Vector denoting current quaternion to base frame as measured by the star tracker.
    Eigen::Vector4d currentQuaternionToBaseFrame_;

    //! Double denoting current altitude as measured by the altimeter.
    double currentAltitude_;

    //! Function to compute the translational acceleration measured by the inertial measurement unit.
    boost::function< void( ) > inertialMeasurementUnitTranslationalAccelerationFunction_;

    //! Function to compute the rotational velocity measured by the inertial measurement unit.
    boost::function< void( ) > inertialMeasurementUnitRotationalVelocityFunction_;

    //! Function to compute the attitude measured by the star tracker.
    boost::function< void( ) > starTrackerOrientationFunction_;

    //! Function to compute the altitude measured by the altimeter.
    boost::function< void( ) > altimeterFunction_;

    //! Function to compute the position and velocity measured by the Deep Space Network.
    boost::function< void( ) > deepSpaceNetworkFunction_;

    //! Map of translational accelerations and rotational velocities measured by the inertial measurment unit
    //! during the current orbit.
    std::map< double, Eigen::Vector6d > currentOrbitHistoryOfInertialMeasurmentUnitMeasurements_;

    //! Map of attitude measured by the star tracker during the current orbit.
    std::map< double, Eigen::Vector4d > currentOrbitHistoryOfStarTrackerMeasurements_;

    //! Map of attitude measured by the star tracker during the current orbit.
    std::map< double, double > currentOrbitHistoryOfAltimeterMeasurements_;

    //! Map of position and velocity measurements by the Deep Space Network during the current simulation.
    std::map< double, std::pair< double, Eigen::Vector6d > > historyOfDeepSpaceNetworkMeasurements_;

    //! Predefined iterator to save (de)allocation time.
    basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationMapIterator_;

    //! Predefined iterator to save (de)allocation time.
    std::map< double, Eigen::Vector6d >::const_iterator inertialMeasurementUnitMeasurementIterator_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_NAVIGATION_INSTRUMENTS_MODEL_H
