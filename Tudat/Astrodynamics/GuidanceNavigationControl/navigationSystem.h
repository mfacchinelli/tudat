/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GNC_NAVIGATION_H
#define TUDAT_GNC_NAVIGATION_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/Filters/filter.h"

namespace tudat
{

namespace guidance_navigation_control
{

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
    NavigationSystem( const boost::shared_ptr< filters::FilterBase< > > navigationFilter,
                      const boost::function< Eigen::MatrixXd( const double, const Eigen::VectorXd& ) >& stateTransitionMatrixFunction,
                      const double planetaryGravitationalParameter,
                      const double planetaryRadius ) :
        navigationFilter_( navigationFilter ), stateTransitionMatrixFunction_( stateTransitionMatrixFunction ),
        planetaryGravitationalParameter_( planetaryGravitationalParameter ), planetaryRadius_( planetaryRadius )
    {
        // Set initial time
        currentTime_ = navigationFilter_->getInitialTime( );

        // Retrieve navigation filter step size
        navigationRefreshStepSize_ = navigationFilter_->getIntegrationStepSize( );

        // Set initial state
        currentEstimatedCartesianState_ = navigationFilter_->getCurrentStateEstimate( );
        currentEstimatedKeplerianState_ =
                orbital_element_conversions::convertCartesianToKeplerianElements( currentEstimatedCartesianState_,
                                                                                  planetaryGravitationalParameter_ );

        // Store initial time and state
        storeCurrentTimeAndStateEstimates( );
    }

    //! Destructor.
    ~NavigationSystem( ) { }

    //! State estimator (SE).
    void stateEstimator( const double previousTime );

    //! Periapse time estimator (PTE).
    /*!
     *  Function to estimate the time of periapsis, by finding the centroid of the aerodynamic acceleration curve,
     *  via the bisection root-finder. Then it computes the change in velocity \f$ \Delta V \f$ by integrating the aerodynamic
     *  acceleration, and the change in semi-major axis, by assuming that the change in velocity occurs as an impulsive shot at
     *  pericenter. Then the values of \f$ \Delta \vartheta \f$ (which is derived from the change in time of periapsis crossing)
     *  and \f$ \Delta a \f$ are used to correct the current state estimate, by using the state transition matrix of the
     *  system.
     *  \param estimatedAerodynamicAcceleration Map of time and estimated aerodynamic acceleration; the acceleration is
     *      the one computed from the IMU measurements and is stored as a three-dimensional vector.
     */
    void periapseTimeEstimator( const std::map< double, Eigen::Vector3d >& estimatedAerodynamicAcceleration );

    //! Atmosphere estimator (AE).
    void atmosphereEstimator( );

    //! Function to retireve current time.
    double getCurrentTime( ) { return currentTime_; }

    //! Function to set current time to new value.
    void setCurrentTime( const double newCurrentTime ) { currentTime_ = newCurrentTime; }

    //! Function to retireve current state.
    std::pair< Eigen::VectorXd, Eigen::VectorXd > getCurrentEstimatedState( )
    {
        return std::make_pair( currentEstimatedCartesianState_, currentEstimatedKeplerianState_ );
    }

    //! Function to set current Cartesian state to new value.
    void setCurrentEstimatedCartesianState( const Eigen::VectorXd& newCurrentCartesianState )
    {
        currentEstimatedCartesianState_ = newCurrentCartesianState;
        currentEstimatedKeplerianState_ =
                orbital_element_conversions::convertCartesianToKeplerianElements( newCurrentCartesianState,
                                                                                  planetaryGravitationalParameter_ );
    }

    //! Function to set current Keplerian state to new value.
    void setCurrentEstimatedKeplerianState( const Eigen::VectorXd& newCurrentKeplerianState )
    {
        currentEstimatedKeplerianState_ = newCurrentKeplerianState;
        currentEstimatedCartesianState_ =
                orbital_element_conversions::convertKeplerianToCartesianElements( newCurrentKeplerianState,
                                                                                  planetaryGravitationalParameter_ );
    }

    //! Function to retrieve the current estimated density.
    double getCurrentEstimatedDensity( ) { return currentEstimatedDensity_; }

    //! Function to retrieve history of estimated states over time.
    /*!
     *  Function to retrieve history of estimated states over time.
     *  \return Map of time and estimated state of the spacecraft; the state is stored as a pair of
     *      Cartesian and Keplerian elements, respectively.
     */
    std::map< double, std::pair< Eigen::VectorXd, Eigen::VectorXd > > getHistoryOfEstimatedStates( )
    {
        return historyOfEstimatedStates_;
    }

    //! Function to retireve standard gravitational parameter of body being orbited.
    double getStandardGravitationalParameter( ) { return planetaryGravitationalParameter_; }

    //! Function to retireve radius of body being orbited.
    double getRadius( ) { return planetaryRadius_; }

private:

    //! Function to store current time and current state estimates.
    void storeCurrentTimeAndStateEstimates( )
    {
        historyOfEstimatedStates_[ currentTime_ ] = std::make_pair( currentEstimatedCartesianState_, currentEstimatedKeplerianState_ );
    }

    //! Double denoting the current time in the estimation process.
    double currentTime_;

    //! Double denoting the integration constant time step for navigation.
    double navigationRefreshStepSize_;

    //! Vector denoting the current estimated state in Cartesian elements.
    Eigen::VectorXd currentEstimatedCartesianState_;

    //! Vector denoting the current estimated state in Keplerian elements.
    Eigen::VectorXd currentEstimatedKeplerianState_;

    //! Double denoting current estimated density.
    double currentEstimatedDensity_;

    //! History of estimated states in Cartesian and Keplerian elements.
    /*!
     *  History of estimated states in Cartesian and Keplerian elements, stored as a map, where the keys are
     *  times and the mapped values are a pair of vectors, denoting the Cartesian and Keplerian elements.
     */
    std::map< double, std::pair< Eigen::VectorXd, Eigen::VectorXd > > historyOfEstimatedStates_;

    //! History of estimated states in Cartesian and Keplerian elements for current orbit.
    /*!
     *  History of estimated states in Cartesian and Keplerian elements for current orbit, stored as a map,
     *  where the keys are times and the mapped values are a pair of vectors, denoting the Cartesian and
     *  Keplerian elements.
     */
    std::map< double, std::pair< Eigen::VectorXd, Eigen::VectorXd > > currentOrbitHistoryOfEstimatedStates_;

    //! Filter object to be used for estimation of state.
    boost::shared_ptr< filters::FilterBase< > > navigationFilter_;

    //! Function to propagate estimated state error.
    boost::function< Eigen::MatrixXd( const double, const Eigen::VectorXd& ) > stateTransitionMatrixFunction_;

    //! Standard gravitational parameter of body being orbited.
    double planetaryGravitationalParameter_;

    //! Radius of body being orbited.
    double planetaryRadius_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_GNC_NAVIGATION_H
