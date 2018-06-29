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

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"

#include "Tudat/Mathematics/Filters/createFilter.h"

//! Typedefs and using statements to simplify code.
typedef Eigen::Matrix< double, 16, 1 > Eigen::Vector16d;
typedef Eigen::Matrix< double, 16, 16 > Eigen::Matrix16d;

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
/*!
 *  Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
 *  The function returns the difference between the areas under the acceleration curve, computed before and after
 *  a true anomaly guess provided by the root-finder object. Since the area under the curve is a constant value, the
 *  slicing procedure will always result in a larger and a smaller value of area, unless the slicing occurs exactly
 *  at the centroid of the area curve (which is the tagert point). Since no interpolation is used, it is possible that
 *  the true anomaly guess is not cointained in the estimatedTrueAnomaly vector, thus, the nearest lower index is used
 *  as reference to compute the area.
 *  \param currentTrueAnomalyGuess Current guess in true anomaly to be used as slicing parameter.
 *  \param estimatedTrueAnomaly Vector of the estimated true anomalies below the atmospheric interface.
 *  \param estimatedAerodynamicAccelerationMagnitude Vector of estimated aerodynamic acceleration magnitudes below the
 *      atmospheric interface altitude.
 *  \return Double representing the difference between the areas under the acceleration curve computed before and
 *      after the current true anomaly estimate (i.e., the slicing parameter).
 */
double areaBisectionFunction( const double currentTrueAnomalyGuess, const std::vector< double >& estimatedTrueAnomaly,
                              const std::vector< double >& estimatedAerodynamicAccelerationMagnitude );

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
                      const double atmosphericInterfaceAltitude ) :
        onboardBodyMap_( onboardBodyMap ), onboardAccelerationModelMap_( onboardAccelerationModelMap ),
        spacecraftName_( spacecraftName ), planetName_( planetName ),
        navigationFilterSettings_( navigationFilterSettings )
    {
        // Set planetary conditions
        planetaryGravitationalParameter_ = onboardBodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( );
        planetaryRadius_ = onboardBodyMap_.at( planetName_ )->getShapeModel( )->getAverageRadius( );
        atmosphericInterfaceRadius_ = planetaryRadius_ + atmosphericInterfaceAltitude;

        // Create root-finder object for bisection of aerodynamic acceleration curve
        // The values inserted are the tolerance in independent value (i.e., one arcminute of true anomaly) and
        // the maximum number of interations (i.e., 10 iterations)
        areaBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 1.7e-3, 10 );
    }

    //! Destructor.
    ~NavigationSystem( ) { }

    //! Function to create navigation filter object for onboard state estimation.
    void createNavigationFilter(
            const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardSystemModel,
            const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardMeasurementModel );

    //! Function to run the State Estimator (SE).
    void runStateEstimator( const double previousTime, const Eigen::Vector4d& currentExternalMeasurementVector );

    //! Function to run the Periapse Time Estimator (PTE).
    /*!
     *  Function to estimate the time of periapsis, by finding the centroid of the aerodynamic acceleration curve,
     *  via the bisection root-finder. Then it computes the change in velocity \f$ \Delta V \f$ by integrating the aerodynamic
     *  acceleration, and the change in semi-major axis, by assuming that the change in velocity occurs as an impulsive shot at
     *  pericenter. Then the values of \f$ \Delta \vartheta \f$ (which is derived from the change in time of periapsis crossing)
     *  and \f$ \Delta a \f$ are used to correct the current state estimate, by using the state transition matrix of the
     *  system.
     *  \param mapOfEstimatedAerodynamicAcceleration Map of time and estimated aerodynamic acceleration; the acceleration is
     *      the one computed from the IMU measurements and is stored as a three-dimensional vector.
     */
    void runPeriapseTimeEstimator( const std::map< double, Eigen::Vector3d >& mapOfEstimatedAerodynamicAcceleration );

    //! Function to run the Atmosphere Estimator (AE).
    void runAtmosphereEstimator( const std::map< double, Eigen::Vector3d >& mapOfEstimatedAerodynamicAcceleration );

    //! Function to update the body and acceleration map with the current time and state information.
    /*!
     *  Function to update the body and acceleration map with the current time and state information.
     *  \param currentTime Time at which the measurement needs to be taken.
     */
    void updateBodyAndAccelerationMaps( const double currentTime )
    {
        // If instruments have not been already updated for the current time
        if ( currentTime_ != currentTime )
        {
            // Update current time
            currentTime_ = currentTime;

            // Update body settings
            onboardBodyMap_.at( spacecraftName_ )->setState( currentEstimatedCartesianState_ );
            onboardBodyMap_.at( spacecraftName_ )->setCurrentRotationalStateToLocalFrame( currentEstimatedRotationalState_ );

            // Update accelerations
            basic_astrodynamics::SingleBodyAccelerationMap accelerationsOnBody = onboardAccelerationModelMap_.at( spacecraftName_ );
            for ( accelerationIterator_ = accelerationsOnBody.begin( ); accelerationIterator_ != accelerationsOnBody.end( );
                  accelerationIterator_++ )
            {
                // Loop over each acceleration
                for ( unsigned int i = 0; i < accelerationIterator_->second.size( ); i++ )
                {
                    // Calculate acceleration and add to state derivative
                    accelerationIterator_->second[ i ]->updateMembers( );
                }
            }
        }
    }

    //! Function to retrieve current estimated translational accelerations exerted on the spacecraft.
    /*!
     *  Function to retrieve current estimated translational accelerations exerted on the spacecraft. The acceleration
     *  is computed by using the
     *  \return
     */
    Eigen::Vector3d getCurrentEstimatedTranslationalAcceleration( )
    {
        // Clear translational accelerations vector
        Eigen::Vector3d currentTranslationalAcceleration;

        // Iterate over all accelerations acting on body
        basic_astrodynamics::SingleBodyAccelerationMap accelerationsOnBody = onboardAccelerationModelMap_.at( spacecraftName_ );
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
        return currentTranslationalAcceleration;
    }

    //! Function to retireve current time.
    double getCurrentTime( ) { return currentTime_; }

    //! Function to retireve current state.
    std::pair< Eigen::Vector6d, Eigen::Vector6d > getCurrentEstimatedTranslationalState( )
    {
        return std::make_pair( currentEstimatedCartesianState_, currentEstimatedKeplerianState_ );
    }

    //! Function to set current Cartesian state to new value.
    void setCurrentEstimatedCartesianState( const Eigen::Vector6d& newCurrentCartesianState )
    {
        currentEstimatedCartesianState_ = newCurrentCartesianState;
        currentEstimatedKeplerianState_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                    newCurrentCartesianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values
    }

    //! Function to set current Keplerian state to new value.
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
     *  \return Map where the keys are times and the mapped value is a pair, whose first element is the pair
     *  of Cartesian and Keplerian elements, whereas the second element is the rotational state.
     */
    std::map< double, std::pair< std::pair< Eigen::Vector6d, Eigen::Vector6d >, Eigen::Vector7d > > getHistoryOfEstimatedStates( )
    {
        return historyOfEstimatedStates_;
    }

    //! Function to retireve the standard gravitational parameter of the body being orbited.
    double getStandardGravitationalParameter( ) { return planetaryGravitationalParameter_; }

    //! Function to retireve the radius of the body being orbited.
    double getRadius( ) { return planetaryRadius_; }

    //! Function to retireve the atmospheric interface radius of the body being orbited.
    double getAtmosphericInterfaceRadius( ) { return atmosphericInterfaceRadius_; }

    //! Clear current orbit estimation history.
    void clearCurrentOrbitEstimationHistory( )
    {
        currentOrbitHistoryOfEstimatedTranslationalStates_.clear( );
    }

private:

    //! Function to store current time and current state estimates.
    void storeCurrentTimeAndStateEstimates( )
    {
        currentOrbitHistoryOfEstimatedTranslationalStates_[ currentTime_ ] = std::make_pair( currentEstimatedCartesianState_,
                                                                                             currentEstimatedKeplerianState_ );
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

    //! Standard gravitational parameter of body being orbited.
    const double planetaryGravitationalParameter_;

    //! Radius of body being orbited.
    const double planetaryRadius_;

    //! Double denoting the atmospheric interface altitude.
    const double atmosphericInterfaceRadius_;

    //! Double denoting the current time in the estimation process.
    double currentTime_;

    //! Double denoting the current time in the estimation process.
    unsigned int currentOrbitCounter_;

    //! Double denoting the integration constant time step for navigation.
    double navigationRefreshStepSize_;

    //! Vector denoting the current estimated translational state in Cartesian elements.
    Eigen::Vector6d currentEstimatedCartesianState_;

    //! Vector denoting the current estimated translational state in Keplerian elements.
    Eigen::Vector6d currentEstimatedKeplerianState_;

    //! Vector denoting the current estimated rotational state.
    Eigen::Vector7d currentEstimatedRotationalState_;

    //! Filter object to be used for estimation of state.
    const boost::shared_ptr< filters::FilterBase< > > navigationFilter_;

    //! Function to propagate estimated state error.
    const boost::function< Eigen::Matrix6d( const Eigen::Vector6d& ) > stateTransitionMatrixFunction_;

    //! Pointer to root-finder used to esimated the time of periapsis.
    boost::shared_ptr< root_finders::BisectionCore< double > > areaBisectionRootFinder_;

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

    //! Predefined iterator to save (de)allocation time.
    basic_astrodynamics::SingleBodyAccelerationMap::const_iterator accelerationIterator_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_NAVIGATION_SYSTEM_H