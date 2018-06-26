/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef MICHELE_GNC_INSTRUMENT
#define MICHELE_GNC_INSTRUMENT

#include <boost/shared_ptr.hpp>
#include <eigen/Core>

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"

namespace tudat
{

namespace guidance_navigation_control
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
class InstrumentSystem
{
public:

    //! Constructor.
    InstrumentSystem( const simulation_setup::NamedBodyMap& bodyMap,
                      const simulation_setup::SelectedAccelerationMap& selectedAccelerationPerBody,
                      const std::map< std::string, std::string >& centralBodies,
                      const std::string& nameOfSpacecraft ) :
        bodyMap_( bodyMap ), selectedAccelerationPerBody_( selectedAccelerationPerBody ),
        centralBodies_( centralBodies ), nameOfSpacecraft_( nameOfSpacecraft )
    {
        // Set instrument presence to false
        inertialMeasurementUnitAdded_ = false;
        starTrackerAdded_ = false;
    }

    //! Destructor.
    ~InstrumentSystem( ) { }

    //! Function to add an inertial measurement unit to the spacecraft set of instruments.
    void addInertialMeasurementUnit( const Eigen::Vector3d& accelerometerBias,
                                     const Eigen::Vector3d& accelerometerScaleFactor,
                                     const Eigen::Vector6d& accelerometerMisalignment,
                                     const Eigen::Vector3d& gyroscopeBias,
                                     const Eigen::Vector3d& gyroscopeScaleFactor,
                                     const Eigen::Vector6d& gyroscopeMisalignment )
    {
        // Check whether an inertial measurement unit is already present
        if ( !inertialMeasurementUnitAdded_ )
        {
            // Inertial measurement unit has been created
            inertialMeasurementUnitAdded_ = true;

            // Create acceleration model object
            accelerationModelMap_ = simulation_setup::createAccelerationModelsMap( bodyMap_, selectedAccelerationPerBody_,
                                                                                   centralBodies_ );

            // Create function for computing corrupted translational accelerations
            inertialMeasurementUnitTranslationalAccelerationFunction_ = boost::bind(
                        &getCurrentTranslationalAcceleration, this, _1, accelerometerBias,
                        computeScaleMisalignmentMatrix( accelerometerScaleFactor, accelerometerMisalignment ) );

            // Create function for computing corrupted rotational velocity
            inertialMeasurementUnitRotationalVelocityFunction_ = boost::bind(
                        &getCurrentRotationalVelocity, this, _1, gyroscopeBias,
                        computeScaleMisalignmentMatrix( gyroscopeScaleFactor, gyroscopeMisalignment ) );
        }
        else
        {
            throw std::runtime_error( "Error in creation of inertial measurement unit for body " + nameOfSpacecraft_ +
                                      ". An IMU is already present." );
        }
    }

    //! Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
    /*!
     *  Function to add a system of orthogonal star trackers to the spacecraft set of instruments.
     *  \param numberOfStarTrackers Number of star trackers to add to the spacecraft.
     *  \param starTrackerAccuracy Accuracy of a star tracker along each axis (3 sigma).
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
            }
            else
            {
                throw std::runtime_error( "Error in creation of star tracker system for body " + nameOfSpacecraft_ +
                                          ". Only a system of two orthogonal star trackers is supported." );
            }
        }
        else
        {
            throw std::runtime_error( "Error in creation of star tracker system for body " + nameOfSpacecraft_ +
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
            Eigen::Vector6d inertialMeasurementUnitMeasurement;

            // Translational accelerations
            Eigen::Vector3d translationalAccelerationVector = Eigen::Vector3d::Zero( );
            inertialMeasurementUnitTranslationalAccelerationFunction_( translationalAccelerationVector );

            // Rotational acceleration (i.e., torque)
            Eigen::Vector3d rotationalVelocityVector = Eigen::Vector3d::Zero( );
            inertialMeasurementUnitRotationalVelocityFunction_( rotationalVelocityVector );

            // Merge translational and rotational accelerations
            inertialMeasurementUnitMeasurement << translationalAccelerationVector, rotationalVelocityVector;
            return inertialMeasurementUnitMeasurement;
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

        }
        else
        {
            throw std::runtime_error( "Error while retrieving rotational velocity from onboard instrument system. No "
                                      "star trackers are present." );
        }
    }

private:

    //! Function to retrieve current translational accelerations exerted on the spacecraft.
    void getCurrentTranslationalAcceleration( Eigen::Vector3d& translationalAcceleration,
                                              const Eigen::Vector3d& biasVector, const Eigen::Matrix3d& scaleMisalignmentMatrix )
    {
        // Iterate over all accelerations acting on body
        basic_astrodynamics::SingleBodyAccelerationMap accelerationsOnBody = accelerationModelMap_.at( nameOfSpacecraft_ );
        for ( accelerationIterator_ = accelerationsOnBody.begin( ); accelerationIterator_ != accelerationsOnBody.end( );
              accelerationIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationIterator_->second.size( ); i++ )
            {
                // Disregard the central gravitational accelerations, since IMUs do not measure them
                if ( i != 0 )
                {
                    // Calculate acceleration and add to state derivative
                    translationalAcceleration += accelerationIterator_->second[ i ]->getAcceleration( );
                }
            }
        }

        // Add errors to acceleration value
        translationalAcceleration = scaleMisalignmentMatrix * translationalAcceleration;
        translationalAcceleration += biasVector;
    }

    //! Function to retrieve current rotational velocity of the spacecraft.
    void getCurrentRotationalVelocity( Eigen::Vector3d& rotationalVelocity,
                                       const Eigen::Vector3d& biasVector, const Eigen::Matrix3d& scaleMisalignmentMatrix )
    {
        // Iterate over all accelerations acting on body
        rotationalVelocity = bodyMap_.at( nameOfSpacecraft_ )->getCurrentAngularVelocityVectorInGlobalFrame( );

        // Add errors to acceleration value
        rotationalVelocity = scaleMisalignmentMatrix * rotationalVelocity;
        rotationalVelocity += biasVector;
    }

    //! Body map of the simulation.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Translational accelerations acting on the spacecraft.
    simulation_setup::SelectedAccelerationMap selectedAccelerationPerBody_;

    //! Central bodies of the simulation.
    std::map< std::string, std::string > centralBodies_;

    //! String denoting the name of the spacecraft body.
    std::string nameOfSpacecraft_;

    //! Boolean denoting whether an inertial measurement unit is present in the spacecraft.
    bool inertialMeasurementUnitAdded_;

    //! Boolean denoting whether a star tracker is present in the spacecraft.
    bool starTrackerAdded_;

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

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // MICHELE_GNC_INSTRUMENT
