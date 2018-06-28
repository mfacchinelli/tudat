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
#include "Tudat/Mathematics/Filters/filter.h"

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
    NavigationSystem( const boost::shared_ptr< filters::FilterBase< > > navigationFilter,
                      const boost::function< Eigen::MatrixXd( const Eigen::VectorXd& ) >& stateTransitionMatrixFunction,
                      const double planetaryGravitationalParameter,
                      const double planetaryRadius,
                      const double atmosphericInterfaceAltitude ) :
        navigationFilter_( navigationFilter ), stateTransitionMatrixFunction_( stateTransitionMatrixFunction ),
        planetaryGravitationalParameter_( planetaryGravitationalParameter ), planetaryRadius_( planetaryRadius ),
        atmosphericInterfaceRadius_( planetaryRadius + atmosphericInterfaceAltitude )
    {
        // Set initial time
        currentTime_ = navigationFilter_->getInitialTime( );
        currentOrbitCounter_ = 0;

        // Retrieve navigation filter step size
        navigationRefreshStepSize_ = navigationFilter_->getIntegrationStepSize( );

        // Set initial state
        currentEstimatedCartesianState_ = navigationFilter_->getCurrentStateEstimate( );
        currentEstimatedKeplerianState_ =
                orbital_element_conversions::convertCartesianToKeplerianElements( currentEstimatedCartesianState_,
                                                                                  planetaryGravitationalParameter_ );

        // Store initial time and state
        storeCurrentTimeAndStateEstimates( );

        // Create root-finder object for bisection of aerodynamic acceleration curve
        // The values inserted are the tolerance in independent value (i.e., one arcminute of true anomaly) and
        //      the maximum number of interations (i.e., 10 iterations)
        areaBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 1.7e-3, 10 );
    }

    //! Destructor.
    ~NavigationSystem( ) { }

    //! State estimator (SE).
    void runStateEstimator( const double previousTime,
                            const boost::shared_ptr< system_models::NavigationInstrumentsModel > instrumentsModel );

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
        currentEstimatedKeplerianState_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                    newCurrentCartesianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values
    }

    //! Function to set current Keplerian state to new value.
    void setCurrentEstimatedKeplerianState( const Eigen::VectorXd& newCurrentKeplerianState )
    {
        currentEstimatedKeplerianState_ = newCurrentKeplerianState;
        currentEstimatedCartesianState_ = orbital_element_conversions::convertKeplerianToCartesianElements(
                    newCurrentKeplerianState, planetaryGravitationalParameter_ );
        storeCurrentTimeAndStateEstimates( ); // overwrite previous values
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

    //! Function to retireve the standard gravitational parameter of the body being orbited.
    double getStandardGravitationalParameter( ) { return planetaryGravitationalParameter_; }

    //! Function to retireve the radius of the body being orbited.
    double getRadius( ) { return planetaryRadius_; }

    //! Function to retireve the atmospheric interface radius of the body being orbited.
    double getAtmosphericInterfaceRadius( ) { return atmosphericInterfaceRadius_; }

    //! Clear current orbit estimation history.
    void clearCurrentOrbitEstimationHistory( )
    {
        currentOrbitHistoryOfEstimatedStates_.clear( );
    }

private:

    //! Function to store current time and current state estimates.
    void storeCurrentTimeAndStateEstimates( )
    {
        historyOfEstimatedStates_[ currentTime_ ] = std::make_pair( currentEstimatedCartesianState_, currentEstimatedKeplerianState_ );
    }

    //! Filter object to be used for estimation of state.
    const boost::shared_ptr< filters::FilterBase< > > navigationFilter_;

    //! Function to propagate estimated state error.
    const boost::function< Eigen::MatrixXd( const Eigen::VectorXd& ) > stateTransitionMatrixFunction_;

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

    //! Vector denoting the current estimated state in Cartesian elements.
    Eigen::VectorXd currentEstimatedCartesianState_;

    //! Vector denoting the current estimated state in Keplerian elements.
    Eigen::VectorXd currentEstimatedKeplerianState_;

    //! Pointer to root-finder used to esimated the time of periapsis.
    boost::shared_ptr< root_finders::BisectionCore< double > > areaBisectionRootFinder_;

    //! History of estimated errors in Keplerian state as computed by the Periapse Time Estimator for each orbit.
    std::map< unsigned int, Eigen::Vector6d > historyOfEstimatedErrorsInKeplerianState_;

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

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_NAVIGATION_SYSTEM_H
