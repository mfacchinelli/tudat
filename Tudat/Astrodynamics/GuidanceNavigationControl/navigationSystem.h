/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NAVIGATION_SYSTEM_H
#define TUDAT_NAVIGATION_SYSTEM_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/customConstantTemperatureAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"

#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"
#include "Tudat/Mathematics/Filters/createFilter.h"

//! Typedefs and using statements to simplify code.
namespace Eigen
{
typedef Eigen::Matrix< double, 16, 1 > Vector16d;
typedef Eigen::Matrix< double, 16, 16 > Matrix16d;
}

namespace tudat
{

namespace guidance_navigation_control
{

//! Enumeration of indices for navigation state vector.
enum NavigationStateIndices
{
    cartesian_position_index = 0,
    cartesian_velocity_index = 3,
    quaternion_real_index = 6,
    quaternion_imaginary_index = 7,
    gyroscope_bias_index = 10,
    gyroscope_scale_factor_index = 13
};

//! Function to remove the error in gyroscope measurement based on the estimated bias and scale factors.
Eigen::Vector3d removeErrorsFromInertialMeasurementUnitMeasurement( const Eigen::Vector3d& currentInertialMeasurementUnitMeasurement,
                                                                    const Eigen::Vector6d& inertialMeasurementUnitErrors );

//! Class for the navigation system of an aerobraking maneuver.
class NavigationSystem
{
public:

    //! Constructor.
    /*!
     *  Constructor for the navigation system of an aerobraking maneuver.
     *  \param navigationFilter
     *  \param stateTransitionMatrixFunction Function to compute the state transition matrix at the current time and state.
     */
    NavigationSystem( const simulation_setup::NamedBodyMap& onboardBodyMap,
                      const basic_astrodynamics::AccelerationMap& onboardAccelerationModelMap,
                      const std::string& spacecraftName, const std::string& planetName,
                      const boost::shared_ptr< filters::FilterSettings< > > navigationFilterSettings,
                      const aerodynamics::AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel,
                      const double atmosphericInterfaceAltitude ) :
        onboardBodyMap_( onboardBodyMap ), onboardAccelerationModelMap_( onboardAccelerationModelMap ),
        spacecraftName_( spacecraftName ), planetName_( planetName ), navigationFilterSettings_( navigationFilterSettings ),
        selectedOnboardAtmosphereModel_( selectedOnboardAtmosphereModel ),
        planetaryGravitationalParameter_( onboardBodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( ) ),
        planetaryRadius_( onboardBodyMap_.at( planetName_ )->getShapeModel( )->getAverageRadius( ) ),
        atmosphericInterfaceRadius_( planetaryRadius_ + atmosphericInterfaceAltitude )
    {
        // Create environment updater
        createOnboardEnvironmentUpdater( );

        // Atmosphere estimator initialization
        atmosphereEstimatorInitialized_ = false;

        // State transition matrix function
        stateTransitionMatrixFunction_ = boost::lambda::constant( Eigen::Matrix6d::Zero( ) );

        // Get index of central body acceleration (which is not measured by the IMUs)
        for ( accelerationMapIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
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
                    break;
                }
            }
        }
    }

    //! Destructor.
    ~NavigationSystem( ) { }

    //! Function to create navigation filter and root-finder objects for onboard state estimation.
    void createNavigationSystemObjects(
            const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardSystemModel,
            const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardMeasurementModel );

    //! Function to run the State Estimator (SE).
    void runStateEstimator( const double currentTime, const Eigen::Vector7d& currentExternalMeasurementVector,
                            const Eigen::Vector3d& currentGyroscopeMeasurement )
    {
        // Set time
        currentTime_ = currentTime;

        // Update filter
        navigationFilter_->updateFilter( currentTime_, currentExternalMeasurementVector );

        // Extract estimated state and update navigation estimates
        Eigen::Vector16d updatedEstimatedState = navigationFilter_->getCurrentStateEstimate( );
        currentEstimatedRotationalState_.segment( 0, 4 ) = updatedEstimatedState.segment( 6, 4 ).normalized( );
        currentEstimatedRotationalState_.segment( 4, 3 ) =
                removeErrorsFromInertialMeasurementUnitMeasurement( currentGyroscopeMeasurement, updatedEstimatedState.segment( 10, 6 ) );
        setCurrentEstimatedCartesianState( updatedEstimatedState.segment( 0, 6 ) );
        // this function also automatically stores the full state estimates at the current time

        // Update body and acceleration maps
        updateOnboardModel( );
    }

    //! Function to run post-atmosphere processes.
    void runPostAtmosphereProcesses( const std::map< double, Eigen::Vector3d >& mapOfEstimatedAerodynamicAcceleration )
    {
        using mathematical_constants::PI;

        // Extract aerodynamic accelerations of when the spacecraft is below the atmospheric interface altitude
        double currentIterationTime;
        std::map< double, Eigen::Vector6d > mapOfEstimatedCartesianStatesBelowAtmosphericInterface;
        std::map< double, Eigen::Vector6d > mapOfEstimatedKeplerianStatesBelowAtmosphericInterface;
        std::map< double, Eigen::Vector3d > mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface;
        std::vector< Eigen::Vector3d > vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface;
        for ( translationalStateIterator_ = currentOrbitHistoryOfEstimatedTranslationalStates_.begin( );
              translationalStateIterator_ != currentOrbitHistoryOfEstimatedTranslationalStates_.end( ); translationalStateIterator_++ )
        {
            if ( translationalStateIterator_->second.first.segment( 0, 3 ).norm( ) <= atmosphericInterfaceRadius_ )
            {
                // Retireve time, state and acceleration of where the altitude is below the atmospheric interface
                currentIterationTime = translationalStateIterator_->first;
                mapOfEstimatedCartesianStatesBelowAtmosphericInterface[ currentIterationTime ] =
                        translationalStateIterator_->second.first;
                mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ] =
                        translationalStateIterator_->second.second;
                mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface[ currentIterationTime ] =
                        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_[ currentIterationTime ];
                vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.push_back(
                            mapOfEstimatedAerodynamicAcceleration.at( currentIterationTime ) );

                // Modify the true anomaly such that it is negative where it is above PI radians (before estimated periapsis)
                if ( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ][ 5 ] >= PI )
                {
                    mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ][ 5 ] -= 2.0 * PI;
                }
            }
        }

        // Remove errors from measured accelerations
        removeAccelerometerErrors( vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface,
                                   mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface );
        // this function also calibrates the accelerometers if it is the first orbit

        // Run post-atmosphere processes if not the first orbit
        if ( true ) // currentOrbitCounter_ != 0 )
        {
            // Compute magnitude of aerodynamic acceleration
            std::vector< double > vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface;
            for ( unsigned int i = 0; i < vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
            {
                vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back(
                            vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ).norm( ) );
            }

            // Run post-atmosphere processes with processed results
            runPeriapseTimeEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                      vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
            runAtmosphereEstimator( mapOfEstimatedCartesianStatesBelowAtmosphericInterface,
                                    vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
        }
    }

    //! Function to retrieve current estimated translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated translational accelerations exerted on the spacecraft. The acceleration
     *  is computed by using the onboard accelerations model map.
     *  \return Vector denoting the full acceleration vector experienced by the spacecraft.
     */
    Eigen::Vector3d getCurrentEstimatedTranslationalAcceleration( )
    {
        // Define output and set accelerations to zero
        Eigen::Vector3d currentTranslationalAcceleration = Eigen::Vector3d::Zero( );

        // Iterate over all accelerations acting on body
        for ( accelerationMapIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each accelerations
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                // Calculate acceleration and add to state derivative
                currentTranslationalAcceleration += accelerationMapIterator_->second[ i ]->getAcceleration( );
            }
        }

        // Add errors to acceleration value
        return currentTranslationalAcceleration;
    }

    //! Function to retrieve current estimated non-gravitational translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated translational accelerations exerted on the spacecraft. The acceleration
     *  is computed by using the onboard accelerations model map. Note that this vector represents only the non-gravitational
     *  accelerations (which are supposed to emulated the accelerations measured by the IMU).
     *  \return Vector denoting only the non-gravitational accelerations.
     */
    Eigen::Vector3d getCurrentEstimatedNonGravitationalTranslationalAcceleration( )
    {
        // Define output and set accelerations to zero
        Eigen::Vector3d currentTranslationalAcceleration = Eigen::Vector3d::Zero( );

        // Iterate over all accelerations acting on body
        for ( accelerationMapIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each accelerations
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                // Disregard the central gravitational accelerations, since IMUs do not measure them
                if ( ( i != sphericalHarmonicsGravityIndex_ ) || ( accelerationMapIterator_->first != planetName_ ) )
                {
                    // Calculate acceleration and add to state derivative
                    currentTranslationalAcceleration += accelerationMapIterator_->second[ i ]->getAcceleration( );
                }
            }
        }

        // Add errors to acceleration value
        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_[ currentTime_ ] = currentTranslationalAcceleration;
        return currentTranslationalAcceleration;
    }

    //! Function to retireve current time.
    double getCurrentTime( ) { return currentTime_; }

    //! Function to retireve current state.
    Eigen::Vector16d getCurrentEstimatedState( )
    {
        return navigationFilter_->getCurrentStateEstimate( );
    }

    //! Function to retireve current state.
    std::pair< Eigen::Vector6d, Eigen::Vector6d > getCurrentEstimatedTranslationalState( )
    {
        return std::make_pair( currentEstimatedCartesianState_, currentEstimatedKeplerianState_ );
    }

    //! Function to compute the current estimated mean motion.
    double getCurrentEstimatedMeanMotion( )
    {
        return basic_astrodynamics::computeKeplerMeanMotion( currentEstimatedKeplerianState_[ 0 ],
                planetaryGravitationalParameter_ );
    }

    //! Function to compute the current estimated mean motion.
    Eigen::Vector6d getEstimatedAccelerometerErrors( )
    {
        return estimatedAccelerometerErrors_;
    }

    //! Function to set current Cartesian state to new value.
    /*!
     *  Function to set current Cartesian state to new value. The value of the Keplerian state is set
     *  automatically by converting the Cartesian state to Keplerian elements. Also, the function updates the
     *  history of translational (and rotational, although unchanged) to the new values.
     *  \param newCurrentKeplerianState New Cartesian state at current time.
     */
    void setCurrentEstimatedCartesianState( const Eigen::Vector6d& newCurrentCartesianState )
    {
        currentEstimatedCartesianState_ = newCurrentCartesianState;
        currentEstimatedKeplerianState_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                    newCurrentCartesianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values
    }

    //! Function to set current Keplerian state to new value.
    /*!
     *  Function to set current Keplerian state to new value. The value of the Cartesian state is set
     *  automatically by converting the Keplerian state to Cartesian elements. Also, the function updates the
     *  history of translational (and rotational, although unchanged) to the new values.
     *  \param newCurrentKeplerianState New Keplerian state at current time.
     */
    void setCurrentEstimatedKeplerianState( const Eigen::Vector6d& newCurrentKeplerianState )
    {
        currentEstimatedKeplerianState_ = newCurrentKeplerianState;
        currentEstimatedCartesianState_ = orbital_element_conversions::convertKeplerianToCartesianElements(
                    newCurrentKeplerianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values
    }

    //! Function to retrieve history of estimated translational and rotational states over time.
    /*!
     *  Function to retrieve history of estimated translational and rotational states over time.
     *  \return Map where the keys are times and the mapped value is a pair, whose first element is the pair of
     *      Cartesian and Keplerian elements, whereas the second element is the rotational state.
     */
    std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > > getHistoryOfEstimatedStates( )
    {
        return historyOfEstimatedStates_;
    }

    //! Function to retrieve history of estimated non-gravitational translational accelerations for the current orbit.
    std::map< double, Eigen::Vector3d > getCurrentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations( )
    {
        return currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_;
    }

    //! Function to retireve refresh step size of navigation system.
    double getNavigationRefreshStepSize( ) { return navigationRefreshStepSize_; }

    //! Function to retireve the standard gravitational parameter of the body being orbited.
    double getStandardGravitationalParameter( ) { return planetaryGravitationalParameter_; }

    //! Function to retireve the radius of the body being orbited.
    double getRadius( ) { return planetaryRadius_; }

    //! Function to retireve the atmospheric interface radius of the body being orbited.
    double getAtmosphericInterfaceRadius( ) { return atmosphericInterfaceRadius_; }

    //! Clear history of estimated states and accelerations for the current orbit.
    void clearCurrentOrbitEstimationHistory( )
    {
        currentOrbitHistoryOfEstimatedTranslationalStates_.clear( );
        currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_.clear( );
    }

    //! Integer denoting the current orbit counter.
    unsigned int currentOrbitCounter_;

private:

    //! Function to create the onboard environment updater.
    void createOnboardEnvironmentUpdater( );

    //! Function to update the body and acceleration map with the current time and state information.
    /*!
     *  Function to update the body and acceleration map with the current time and state information.
     */
    void updateOnboardModel( )
    {
        // Update environment
        std::unordered_map< propagators::IntegratedStateType, Eigen::VectorXd > mapOfStatesToUpdate;
        mapOfStatesToUpdate[ propagators::translational_state ] = currentEstimatedCartesianState_ +
                onboardBodyMap_.at( planetName_ )->getState( );
        mapOfStatesToUpdate[ propagators::rotational_state ] = currentEstimatedRotationalState_;
        onboardEnvironmentUpdater_->updateEnvironment( currentTime_, mapOfStatesToUpdate );

//        // Update body settings
//        onboardBodyMap_.at( spacecraftName_ )->setState( currentEstimatedCartesianState_ );
//        onboardBodyMap_.at( spacecraftName_ )->setCurrentRotationalStateToLocalFrame( currentEstimatedRotationalState_ );

        // Loop over bodies exerting accelerations on spacecraft
        for ( accelerationMapIterator_ = onboardAccelerationModelMap_.at( spacecraftName_ ).begin( );
              accelerationMapIterator_ != onboardAccelerationModelMap_.at( spacecraftName_ ).end( );
              accelerationMapIterator_++ )
        {
            // Loop over each acceleration
            for ( unsigned int i = 0; i < accelerationMapIterator_->second.size( ); i++ )
            {
                // Update acceleration model
                accelerationMapIterator_->second[ i ]->updateMembers( currentTime_ );
            }
        }
    }

    //! Function to remove and calibrate (first time only) accelerometer errors.
    /*!
     *  Function to remove and calibrate (first time only) accelerometer errors.
     *  \param vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface Vector of estimated aerodynamic acceleration
     *      below the atmospheric interface altitude and corrupted by accelerometer errors. The acceleration is thus taken
     *      directly from the IMU.
     *  \param mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface Map of time and estimated Cartesian elements below
     *      the atmospheric interface altitude.
     */
    void removeAccelerometerErrors(
            std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface,
            const std::map< double, Eigen::Vector3d >& mapOfExpectedAerodynamicAccelerationBelowAtmosphericInterface );

    //! Function to run the Periapse Time Estimator (PTE).
    /*!
     *  Function to estimate the time of periapsis, by finding the centroid of the aerodynamic acceleration curve,
     *  via the bisection root-finder. Then it computes the change in velocity \f$ \Delta V \f$ by integrating the aerodynamic
     *  acceleration, and the change in semi-major axis, by assuming that the change in velocity occurs as an impulsive shot at
     *  pericenter. Then the values of \f$ \Delta \vartheta \f$ (which is derived from the change in time of periapsis crossing)
     *  and \f$ \Delta a \f$ are used to correct the current state estimate, by using the state transition matrix of the
     *  system.
     *  \param mapOfEstimatedKeplerianStatesBelowAtmosphericInterface Map of time and estimated Keplerian elements below
     *      the atmospheric interface altitude.
     *  \param mapOfEstimatedAerodynamicAcceleration Map of time and estimated aerodynamic acceleration. The acceleration is
     *      the one computed from the IMU measurements and is stored as a three-dimensional vector.
     */
    void runPeriapseTimeEstimator(
            const std::map< double, Eigen::Vector6d >& mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
            const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

    //! Function to run the Atmosphere Estimator (AE).
    /*!
     *  Function to estimate the atmospheric properties based on the selected onboard atmosphere model.
     *  \param mapOfEstimatedCartesianStatesBelowAtmosphericInterface Map of time and estimated Cartesian elements below
     *      the atmospheric interface altitude.
     *  \param mapOfEstimatedAerodynamicAcceleration Map of time and estimated aerodynamic acceleration. The acceleration is
     *      the one computed from the IMU measurements and is stored as a three-dimensional vector.
     */
    void runAtmosphereEstimator(
            const std::map< double, Eigen::Vector6d >& mapOfEstimatedCartesianStatesBelowAtmosphericInterface,
            const std::vector< double >& vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

    //! Function to store current time and current state estimates.
    void storeCurrentTimeAndStateEstimates( )
    {
        // Store translational values for current orbit
        currentOrbitHistoryOfEstimatedTranslationalStates_[ currentTime_ ] = std::make_pair( currentEstimatedCartesianState_,
                                                                                             currentEstimatedKeplerianState_ );

        // Store full state for navigation history
        historyOfEstimatedStates_[ currentTime_ ] = std::make_pair( std::make_pair( currentEstimatedCartesianState_,
                                                                                    currentEstimatedKeplerianState_ ),
                                                                    currentEstimatedRotationalState_ );
    }

    //! Body map of the simulation.
    const simulation_setup::NamedBodyMap onboardBodyMap_;

    //! Pointer to accelerations exerted on the spacecraft.
    const basic_astrodynamics::AccelerationMap onboardAccelerationModelMap_;

    //! String denoting the name of the spacecraft body.
    const std::string spacecraftName_;

    //! String denoting the name of the planet being orbited body.
    const std::string planetName_;

    //! Pointer to the filter settings to be used to create the navigation filter for state estimation.
    const boost::shared_ptr< filters::FilterSettings< > > navigationFilterSettings_;

    //! Enumeration denoting atmosphere model being used for onboard applications.
    const aerodynamics::AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel_;

    //! Standard gravitational parameter of body being orbited.
    const double planetaryGravitationalParameter_;

    //! Radius of body being orbited.
    const double planetaryRadius_;

    //! Double denoting the atmospheric interface altitude.
    const double atmosphericInterfaceRadius_;

    //! Pointer to the onboard environment updater.
    /*!
     *  Pointer to the onboard environment updater. The environment is updated based on the current state and time.
     *  Calling the updateEnvironment function automatically updates all dependent variables that are needed to calulate
     *  the state derivative.
     */
    boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > onboardEnvironmentUpdater_;

    //! Double denoting the current time in the estimation process.
    double currentTime_;

    //! Double denoting the integration constant time step for navigation.
    double navigationRefreshStepSize_;

    //! Vector denoting the current estimated translational state in Cartesian elements.
    Eigen::Vector6d currentEstimatedCartesianState_;

    //! Vector denoting the current estimated translational state in Keplerian elements.
    Eigen::Vector6d currentEstimatedKeplerianState_;

    //! Vector denoting the current estimated rotational state.
    /*!
     *  Vector denoting the current estimated rotational state, where the quaternion expresses the rotation from
     *  body-fixed (local) to inertial (global or base) frame, whereas the rotational velocity is expressed in the body frame.
     */
    Eigen::Vector7d currentEstimatedRotationalState_;

    //! Filter object to be used for estimation of state.
    boost::shared_ptr< filters::FilterBase< > > navigationFilter_;

    //! Pointer to root-finder used to esimated the time of periapsis.
    boost::shared_ptr< root_finders::BisectionCore< double > > areaBisectionRootFinder_;

    //! Function to propagate estimated state error.
    boost::function< Eigen::Matrix6d( const Eigen::Vector6d& ) > stateTransitionMatrixFunction_;

    //! Integer denoting the index of the spherical harmonics gravity exerted by the planet being orbited.
    unsigned int sphericalHarmonicsGravityIndex_;

    //! Vector denoting the bias and scale errors of the accelerometer after calibration.
    Eigen::Vector6d estimatedAccelerometerErrors_;

    //! Boolean denoting whether the atmosphere estimator has been initialized.
    /*!
     *  Boolean denoting whether the atmosphere estimator has been initialized. Note that the initialization is assumed to
     *  be achieved once at least 7 orbits have been carried out.
     */
    bool atmosphereEstimatorInitialized_;

    //! History of estimated errors in Keplerian state as computed by the Periapse Time Estimator for each orbit.
    std::map< unsigned int, Eigen::Vector6d > historyOfEstimatedErrorsInKeplerianState_;

    //! History of estimated translational (in Cartesian and Keplerian elements) and rotational states.
    /*!
     *  History of estimated translational (in Cartesian and Keplerian elements) and rotational states,
     *  stored as a map, where the keys are times and the mapped value is a pair, whose first element is the pair
     *  of Cartesian and Keplerian elements, whereas the second element is the rotational state.
     */
    std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > > historyOfEstimatedStates_;

    //! History of estimated states in Cartesian and Keplerian elements for current orbit.
    /*!
     *  History of estimated states in Cartesian and Keplerian elements for current orbit, stored as a map,
     *  where the keys are times and the mapped values are a pair of vectors, denoting the Cartesian and
     *  Keplerian elements.
     */
    std::map< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > > currentOrbitHistoryOfEstimatedTranslationalStates_;

    //! History of estimated non-gravitational translational accelerations for current orbit.
    std::map< double, Eigen::Vector3d > currentOrbitHistoryOfEstimatedNonGravitationalTranslationalAccelerations_;

    //! Predefined iterator to save (de)allocation time.
    basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationMapIterator_;

    //! Predefined iterator to save (de)allocation time.
    std::map< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > >::const_iterator translationalStateIterator_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_NAVIGATION_SYSTEM_H
