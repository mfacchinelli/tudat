/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GNC_EXTRA_FUNCTIONS_H
#define TUDAT_GNC_EXTRA_FUNCTIONS_H

#include <map>
#include <tuple>
#include <vector>
#include <iostream>

#include <boost/function.hpp>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"
#include "Tudat/Mathematics/RootFinders/bisection.h"

namespace tudat
{

namespace system_models
{

//! Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
/*!
 *  Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
 *  The function returns the difference between the areas under the acceleration curve, computed before and after
 *  a time guess provided by the root-finder object. Since the area under the curve is a constant value, the
 *  slicing procedure will always result in a larger and a smaller value of area, unless the slicing occurs exactly
 *  at the centroid of the area curve (which is the tagert point). Since no interpolation is used, it is possible that
 *  the time guess is not cointained in the onboardTime vector, thus, the nearest index is used as reference to
 *  compute the area.
 *  \param currentTimeGuess Current guess in periapse time to be used as slicing parameter.
 *  \param constantTimeStep Value of constant time step to be used for numerical integration.
 *  \param onboardTime Onboard times below the atmospheric interface.
 *  \param estimatedAerodynamicAccelerationMagnitude Vector of estimated aerodynamic acceleration magnitudes below the
 *      atmospheric interface altitude.
 *  \return Double representing the difference between the areas under the acceleration curve computed before and
 *      after the current time estimate (i.e., the slicing parameter).
 */
double areaBisectionFunction( const double currentTimeGuess, const double constantTimeStep,
                              const Eigen::VectorXd& onboardTime,
                              const std::vector< double >& estimatedAerodynamicAccelerationMagnitude );

//! Function to be used as input to the root-finder to determine the magnitude of the apoapsis maneuver.
/*!
 *  Function to be used as input to the root-finder to determine the magnitude of the apoapsis maneuver.
 *  \param currentMagnitudeGuess
 *  \param initialEstimatedCartesianState
 *  \param targetPeriapsisRadius
 *  \param transformationFromLocalToInertialFrame
 *  \param statePropagationFunction
 *  \return
 */
double maneuverBisectionFunction( const double currentMagnitudeGuess, const Eigen::Vector6d& initialEstimatedCartesianState,
                                  const double targetPeriapsisRadius, const Eigen::Matrix3d& transformationFromLocalToInertialFrame,
                                  const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
                                  std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction );

//! Function to be used as input to the non-linear least squares process to determine the accelerometer errors.
/*!
 *  Function to be used as input to the non-linear least squares process to determine the accelerometer errors.
 *  \param currentErrorEstimate
 *  \param vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface Vector of estimated acceleration
 *  \return Pair of expected acceleration and Jacobian of the acceleration function w.r.t. the accelerometer errors.
 */
std::pair< Eigen::VectorXd, Eigen::MatrixXd > accelerometerErrorEstimationFunction(
        const Eigen::Vector6d& currentErrorEstimate,
        const std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface );

//! Function to be used as input to the non-linear least squares process to determine the atmosphere parameters of the
//! three-term atmosphere model.
/*!
 *  Function to be used as input to the non-linear least squares process to determine the atmosphere parameters of the
 *  three-term atmosphere model.
 *  \param currentParameterEstimate
 *  \param vectorOfEstimatedAltitudesBelowAtmosphericInterface
 *  \param referenceAltitude
 *  \return
 */
std::pair< Eigen::VectorXd, Eigen::MatrixXd > threeModelParametersEstimationFunction(
        const Eigen::Vector5d& currentParameterEstimate,
        const Eigen::VectorXd& vectorOfEstimatedAltitudesBelowAtmosphericInterface,
        const double referenceAltitude );

//! Altitude correction function for guidance system corridor estimator.
double correctionFactorForCorridorBoundaries( const double altitudeGuess, const Eigen::Vector2d& linearLeastSquaresEstimate );

//! Class for corridor estimator of an aerobraking maneuver.
class CorridorEstimator
{
public:

    //! Enumeration of possible boundary types.
    enum BoundaryType
    {
        lower_heating = 0,
        lower_lifetime = 1,
        upper = 2
    };

    //! Constructor.
    CorridorEstimator( const double maximumHeatRate, const double maximumHeatLoad,
                       const double minimumDynamicPressure, const double minimumLifetime,
                       const std::pair< double, double >& boundariesForLowerAltitudeBasedOnHeating,
                       const std::pair< double, double >& boundariesForLowerAltitudeBasedOnLifetime,
                       const std::pair< double, double >& boundariesForUpperAltitude,
                       const double planetaryGravitationalParameter, const double planetaryRadius ) :
        maximumHeatRate_( maximumHeatRate ), maximumHeatLoad_( maximumHeatLoad ),
        minimumDynamicPressure_( minimumDynamicPressure ), minimumLifetime_( minimumLifetime ),
        boundariesForLowerAltitudeBasedOnHeating_( boundariesForLowerAltitudeBasedOnHeating ),
        boundariesForLowerAltitudeBasedOnLifetime_( boundariesForLowerAltitudeBasedOnLifetime ),
        boundariesForUpperAltitude_( boundariesForUpperAltitude ),
        planetaryGravitationalParameter_( planetaryGravitationalParameter ), planetaryRadius_( planetaryRadius )
    {
        // Create root-finder object for bisection of periapsis altitude
        // The values inserted are the tolerance in independent value (i.e., the percentage corresponding to 0.5 km difference at
        // 100 km altitude) and the maximum number of iterations (i.e., 10 iterations)
        altitudeBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 0.5 / 100.0, 10 );
    }

    //! Destructor.
    ~CorridorEstimator( ) { }

    //! Function to run the estimator for the specified corridor boundary.
    double estimateCorridorBoundary(
            const BoundaryType typeOfBoundaryToEstimate, const Eigen::Vector6d& initialEstimatedKeplerianState,
            const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
            std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction );

    //! Function to retireve the function returning the multiplication factor to be used to correct for the difference between the
    //! Kepler orbit assumption of the corridor estimator.
    /*!
     *  Function to retireve the function returning the multiplication factor to be used to correct for the difference between the
     *  Kepler orbit assumption of the corridor estimator.
     *  \return Double-returning function, where the output denotes the multiplication factor for the corridor boundaries and the input
     *      represents the estimated corridor boundary value.
     */
    boost::function< double( const double ) > getPeriapsisAltitudeCorrectionFunction( )
    {
        // Extract altitude information
        Eigen::VectorXd vectorOfAltitudeGuesses;
        vectorOfAltitudeGuesses.resize( historyOfPeriapsisInformation_.size( ) );
        Eigen::VectorXd vectorOfAltitudeRatios;
        vectorOfAltitudeRatios.resize( historyOfPeriapsisInformation_.size( ) );
        for ( unsigned int i = 0; i < historyOfPeriapsisInformation_.size( ); i++ )
        {
            vectorOfAltitudeGuesses[ i ] = historyOfPeriapsisInformation_.at( i ).first;
            vectorOfAltitudeRatios[ i ] = vectorOfAltitudeGuesses[ i ] / historyOfPeriapsisInformation_.at( i ).second;
        }

        // Clear history for next estimation
        historyOfPeriapsisInformation_.clear( );

        // Use least squares to estimate value of linear regression
        std::vector< double > vectorOfPolynomialPowers = { 0, 1 };
        Eigen::Vector2d estimatedLinearCoefficients = linear_algebra::getLeastSquaresPolynomialFit(
                    vectorOfAltitudeGuesses, vectorOfAltitudeRatios, vectorOfPolynomialPowers );

        // Give output
        return boost::bind( &correctionFactorForCorridorBoundaries, _1, estimatedLinearCoefficients );
    }

private:

    //! Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
    /*!
     *  Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor. This function
     *  uses the maximum allowed heat rate and load as constraints.
     *  \param currentAltitudeGuess Double denoting the current altitude guess of the bisection root-finder.
     *  \param initialEstimatedKeplerianState Vector denoting the initial value of estimated Keplerian elements.
     *  \param statePropagationFunction Function used to propagate the spacecraft position, based on custom termination settings
     *      and custom initial conditions. The output of the function is a pair, where the first element is a boolean denoting
     *      whether the propagation was successful, and the second element is another pair, where the first entry is the state
     *      history and the second entry the dependent variable history.
     *  \return Double denoting the value of the lower altitude bound for the periapsis corridor, estimated by using heating conditions
     *      as active constraints.
     */
    double lowerAltitudeBisectionFunctionBasedOnHeatingConditions(
            const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
            const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
            std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction );

    //! Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
    /*!
     *  Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor. This function
     *  uses the minimum allowed lifetime as constraint.
     *  \param currentAltitudeGuess Double denoting the current altitude guess of the bisection root-finder.
     *  \param initialEstimatedKeplerianState Vector denoting the initial value of estimated Keplerian elements.
     *  \param statePropagationFunction Function used to propagate the spacecraft position, based on custom termination settings
     *      and custom initial conditions. The output of the function is a pair, where the first element is a boolean denoting
     *      whether the propagation was successful, and the second element is another pair, where the first entry is the state
     *      history and the second entry the dependent variable history.
     *  \return Double denoting the value of the lower altitude bound for the periapsis corridor, estimated by using lifetime as active
     *      constraint.
     *  \return
     */
    double lowerAltitudeBisectionFunctionBasedOnLifetimeCondition(
            const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
            const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
            std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction );

    //! Function to be used as input to the root-finder to determine the upper altitude bound for the periapsis corridor.
    /*!
     *  Function to be used as input to the root-finder to determine the upper altitude bound for the periapsis corridor.
     *  \param currentAltitudeGuess Double denoting the current altitude guess of the bisection root-finder.
     *  \param initialEstimatedKeplerianState Vector denoting the initial value of estimated Keplerian elements.
     *  \param statePropagationFunction Function used to propagate the spacecraft position, based on custom termination settings
     *      and custom initial conditions. The output of the function is a pair, where the first element is a boolean denoting
     *      whether the propagation was successful, and the second element is another pair, where the first entry is the state
     *      history and the second entry the dependent variable history.
     *  \return Double denoting the value of the upper altitude bound for the periapsis corridor.
     *  \return
     */
    double upperAltitudeBisectionFunction( const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
                                           const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
                                           std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction );

    //! Function to propagate the state based on the current altitude guess for periapsis as initial condition.
    /*!
     *  \brief propagateStateWithAltitudeGuess
     *  \param currentAltitudeGuess Double denoting the current altitude guess of the bisection root-finder.
     *  \param initialEstimatedKeplerianState Vector denoting the initial value of estimated Keplerian elements.
     *  \param statePropagationFunction Function used to propagate the spacecraft position, based on custom termination settings
     *      and custom initial conditions. The output of the function is a pair, where the first element is a boolean denoting
     *      whether the propagation was successful, and the second element is another pair, where the first entry is the state
     *      history and the second entry the dependent variable history.
     *  \return Double denoting the value of the lower altitude bound for the periapsis corridor.
     *  \return Pair of propagation history, where the first entry is the state history and the second entry the dependent variable history.
     */
    std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagateStateWithAltitudeGuess(
            const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
            const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
            std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
    {
        // Copy initial Keplerian state
        Eigen::Vector6d initialKeplerianState = initialEstimatedKeplerianState;

        // Modify initial state to match new estimated periapsis altitude
        double estimatedApoapsisRadius = basic_astrodynamics::computeKeplerRadialDistance(
                    initialKeplerianState[ 0 ], initialKeplerianState[ 1 ], initialKeplerianState[ 5 ] );
        double semiMajorAxis = 0.5 * ( estimatedApoapsisRadius + currentAltitudeGuess + planetaryRadius_ );
        double eccentricity = ( estimatedApoapsisRadius - currentAltitudeGuess - planetaryRadius_ ) /
                ( estimatedApoapsisRadius + currentAltitudeGuess + planetaryRadius_ );
        initialKeplerianState[ 0 ] = semiMajorAxis;
        initialKeplerianState[ 1 ] = eccentricity;

        // Propagate orbit to new condition and retrieve heating conditions
        return statePropagationFunction( orbital_element_conversions::convertKeplerianToCartesianElements(
                                              initialKeplerianState, planetaryGravitationalParameter_ ) ).second;
    }

    //! Double denoting the maximum allowed heat flux that the spacecraft can endure.
    const double maximumHeatRate_;

    //! Double denoting the maximum allowed heat load that the spacecraft can endure.
    const double maximumHeatLoad_;

    //! Double denoting the miminum allowed dynamic pressure that the spacecraft should encounter.
    const double minimumDynamicPressure_;

    //! Double denoting the minimum allowed predicted lifetime in days.
    const double minimumLifetime_;

    //! Pair denoting the lower and upper bound for the lower altitude root finder based on heating conditions.
    const std::pair< double, double > boundariesForLowerAltitudeBasedOnHeating_;

    //! Pair denoting the lower and upper bound for the lower altitude root finder based on lifetime conditions.
    const std::pair< double, double > boundariesForLowerAltitudeBasedOnLifetime_;

    //! Pair denoting the lower and upper bound for the upper altitude root finder.
    const std::pair< double, double > boundariesForUpperAltitude_;

    //! Standard gravitational parameter of body being orbited.
    const double planetaryGravitationalParameter_;

    //! Radius of body being orbited.
    const double planetaryRadius_;

    //! Pointer to root-finder used to esimate the periapsis corridor altitudes.
    boost::shared_ptr< root_finders::BisectionCore< double > > altitudeBisectionRootFinder_;

    //! Vector of pairs denoting the root-finder altitude guess and the actual (propagated) periapsis altitude.
    std::vector< std::pair< double, double > > historyOfPeriapsisInformation_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_GNC_EXTRA_FUNCTIONS_H
