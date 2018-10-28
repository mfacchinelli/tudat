/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Facchinelli, M. (2018). Aerobraking Navigation, Guidance and Control.
 *          Master Thesis, Delft University of Technology.
 */

#ifndef TUDAT_GUIDANCE_SYSTEM_H
#define TUDAT_GUIDANCE_SYSTEM_H

#include <map>
#include <tuple>
#include <vector>

#include <boost/lambda/lambda.hpp>

#include <Eigen/Geometry>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Basics/utilities.h"

#include "Tudat/Astrodynamics/SystemModels/extraFunctions.h"
#include "Tudat/Mathematics/RootFinders/bisection.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"

namespace tudat
{

namespace system_models
{

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
        // The values inserted are the tolerance in independent value (i.e., the percentage corresponding to 100 m difference at
        // 100 km altitude) and the maximum number of iterations (i.e., 10 iterations)
        altitudeBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 0.1 / 100.0, 10 );
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

//! Class for guidance system of an aerobraking maneuver.
class GuidanceSystem
{
public:

    //! Enumeration for aerobraking phases.
    enum AerobrakingPhaseIndicator
    {
        walk_in_phase = 0,
        main_phase = 1,
        walk_out_phase = 2,
        termination_phase = 3,
        aerobraking_complete = 4
    };

    //! Constructor.
    /*!
     *  Constructor.
     *  \param targetPeriapsisAltitude Double denoting the value of the target periapsis altitude.
     *  \param targetApoapsisAltitude Double denoting the value of the target apoapsis altitude.
     *  \param maximumAllowedHeatRate Double denoting the maximum allowed heat rate that the spacecraft can endure.
     *  \param maximumAllowedHeatLoad Double denoting the maximum allowed heat load that the spacecraft can endure.
     *  \param minimumAllowedDynamicPressure Double denoting the minimum allowed dynamic pressure that the spacecraft should encounter.
     *  \param minimumAllowedLifetime Double denoting the minimum allowed predicted lifetime in days.
     *  \param guidanceTesting Boolean denoting whether the guidance system is being tested.
     */
    GuidanceSystem( const double targetPeriapsisAltitude,
                    const double targetApoapsisAltitude,
                    const double maximumAllowedHeatRate,
                    const double maximumAllowedHeatLoad,
                    const double minimumAllowedDynamicPressure,
                    const double minimumAllowedLifetime,
                    const bool guidanceTesting = false ) :
        targetPeriapsisAltitude_( targetPeriapsisAltitude ), targetApoapsisAltitude_( targetApoapsisAltitude ),
        maximumAllowedHeatRate_( maximumAllowedHeatRate ), maximumAllowedHeatLoad_( maximumAllowedHeatLoad ),
        minimumAllowedDynamicPressure_( minimumAllowedDynamicPressure ), minimumAllowedLifetime_( minimumAllowedLifetime ),
        guidanceTesting_( guidanceTesting )
    {
        // Inform user
        if ( guidanceTesting_ )
        {
            std::cerr << "Warning in guidance system. Guidance system testing is active." << std::endl;
        }

        // Create root-finder object for bisection of maneuver magnitude estimate
        // The values inserted are the tolerance in independent value (i.e., the percentage corresponding to 100 m difference at
        // 100 km altitude) and the maximum number of iterations (i.e., 10 iterations)
        maneuverBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 0.1 / 100.0, 10 );

        // Set values to their initial conditions
        periapsisAltitudeScaling_ = TUDAT_NAN;
    }

    //! Destructor.
    ~GuidanceSystem( ) { }

    //! Function to create the guidance system objects.
    /*!
     *  Function to create the guidance system objects.
     *  \param statePropagationFunction Function used to propagate the spacecraft position, based on custom termination settings
     *      and custom initial conditions. The output of the function is a pair, where the first element is a boolean denoting
     *      whether the propagation was successful, and the second element is another pair, where the first entry is the state
     *      history and the second entry the dependent variable history.
     *  \param planetaryGravitationalParameter Double denoting the gravitational parameter of the planet being orbited.
     *  \param planetaryRadius Double denoting the radius of the planet being orbited.
     */
    void createGuidanceSystemObjects(
            const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
            std::map< double, Eigen::VectorXd > > >( const boost::shared_ptr< propagators::PropagationTerminationSettings >,
                const Eigen::Vector6d& ) >& statePropagationFunction,
            const double planetaryGravitationalParameter, const double planetaryRadius )
    {
        // Set planet parameters
        planetaryGravitationalParameter_ = planetaryGravitationalParameter;
        planetaryRadius_ = planetaryRadius;

        // Create propagation function
        statePropagationFunction_ = statePropagationFunction;

        // Create corridor estimator
        corridorEstimator_ = boost::make_shared< CorridorEstimator >( maximumAllowedHeatRate_, maximumAllowedHeatLoad_,
                                                                      minimumAllowedDynamicPressure_, minimumAllowedLifetime_,
                                                                      std::make_pair( 90.0e3, 130.0e3 ),
                                                                      std::make_pair( 100.0e3, 140.0e3 ),
                                                                      std::make_pair( 100.0e3, 130.0e3 ),
                                                                      planetaryGravitationalParameter_, planetaryRadius_ );
    }

    //! Function to determine in which aerobraking phase the spacecraft is currently in.
    /*!
     *  Function to determine in which aerobraking phase the spacecraft is currently in. Aerobraking is divided in four phases:
     *  walk-in, main, walk-out and completed. During the walk-in phase, the periapsis is slowly lowered into the atmosphere,
     *  \param currentEstimatedKeplerianState Vector denoting the current estimated translational Keplerian elements.
     *  \param pairOfAtmosphereInitiationIndicators Pair of integers, where the first element denotes the number of atmosphere
     *      samples that have been taken so far, and the second element indicates the number of atmosphere samples required for
     *      the atmosphere estimator to be considered initialized.
     */
    void determineAerobrakingPhase( const Eigen::Vector6d& currentEstimatedKeplerianState,
                                    const std::pair< unsigned int, unsigned int >& pairOfAtmosphereInitiationIndicators )
    {
        // Store previous aerobraking phase indicator
        AerobrakingPhaseIndicator previousAerobrakingPhase = currentOrbitAerobrakingPhase_;

        // Declare current orbit aerobraking phase indicator and set value to main phase
        AerobrakingPhaseIndicator detectedAerobrakingPhase = main_phase;

        // Set periapsis altitude scaling to default value
        periapsisAltitudeScaling_ = 1.0;

        // Check whether atmosphere has been initiated
        if ( pairOfAtmosphereInitiationIndicators.first < pairOfAtmosphereInitiationIndicators.second )
        {
            detectedAerobrakingPhase = walk_in_phase;
            periapsisAltitudeScaling_ = 1.2 - 0.2 * ( static_cast< double >( pairOfAtmosphereInitiationIndicators.first ) /
                                                      static_cast< double >( pairOfAtmosphereInitiationIndicators.second ) );
        }

        // Check whether apoapsis is approaching target value
        if ( currentEstimatedKeplerianState[ 1 ] < 0.3 )
        {
            detectedAerobrakingPhase = walk_out_phase;
        }

        // Check whether it is time to perform periapsis raise maneuver is complete
        if ( ( computeCurrentFirstOrderEstimatedApoapsisRadius( currentEstimatedKeplerianState ) -
               planetaryRadius_ ) < 1.25 * targetApoapsisAltitude_ )
        {
            detectedAerobrakingPhase = termination_phase;
        }

        // Check whether aerobraking is complete
        if ( previousAerobrakingPhase == termination_phase )
        {
            detectedAerobrakingPhase = aerobraking_complete;
        }

        // Set current orbit aerobraking phase
        currentOrbitAerobrakingPhase_ = detectedAerobrakingPhase;
    }

    //! Function to run corridor estimator (CE).
    /*!
     *  Function to run corridor estimator (CE). This function estimates the lower and upper altitude bounds for the periapsis
     *  corridor. The lower bound corresponds to the altitude where the estimated heat rate and heat load are below the maximum
     *  allowed value (which depend on the spacecraft material properties), whereas the upper bound corresponds to the altitude
     *  where the dynamic pressure is higher than the minimum allowed value (which depends on the total aerobraking duration).
     *  \param currentTime Double denoting the current time.
     *  \param currentEstimatedCartesianState Vector denoting the current estimated translational Cartesian elements.
     *  \param currentEstimatedKeplerianState Vector denoting the current estimated translational Keplerian elements.
     */
    void runCorridorEstimator( const double currentTime,
                               const Eigen::Vector6d& currentEstimatedCartesianState,
                               const Eigen::Vector6d& currentEstimatedKeplerianState );

    //! Function to run apoapsis maneuver estimator (ME).
    /*!
     *  Function to run apoapsis maneuver estimator (ME). Note that since the root-finder uses the propagator defined in the
     *  runCorridorEstimator function (i.e., periodReducedStatePropagationFunction_), said function needs to be run before this one.
     *  \param currentEstimatedCartesianState Vector denoting the current estimated translational Cartesian elements.
     *  \param currentEstimatedKeplerianState Vector denoting the current estimated translational Keplerian elements.
     *  \param currentEstimatedMeanMotion Double denoting the current estimated mean motion.
     *  \param improveEstimateWithBisection Boolean denoting whether the maneuver estimate should be improved by using a
     *      bisection root-finder algorithm.
     */
    void runApoapsisManeuverEstimator( const Eigen::Vector6d& currentEstimatedCartesianState,
                                       const Eigen::Vector6d& currentEstimatedKeplerianState,
                                       const double currentEstimatedMeanMotion,
                                       const bool improveEstimateWithBisection = true );

    //! Function to run periapsis maneuver estimator (ME).
    /*!
     *  Function to run periapsis maneuver estimator (ME).
     *  \param currentTime Double denoting the current time.
     *  \param currentEstimatedCartesianState Vector denoting the current estimated translational Cartesian elements.
     *  \param currentEstimatedKeplerianState Vector denoting the current estimated translational Keplerian elements.
     *  \param currentEstimatedMeanMotion Double denoting the current estimated mean motion.
     */
    void runPeriapsisManeuverEstimator( const double currentTime,
                                        const Eigen::Vector6d& currentEstimatedCartesianState,
                                        const Eigen::Vector6d& currentEstimatedKeplerianState,
                                        const double currentEstimatedMeanMotion );

    //! Function to set the value of the current orbit counter.
    void setCurrentOrbitCounter( const unsigned int currentOrbitCounter )
    {
        currentOrbitCounter_ = currentOrbitCounter;
    }

    //! Function to retrieve whether the apoapsis maneuver is to be performed.
    /*!
     *  Function to retrieve whether the apoapsis maneuver is to be performed. Note that the inverse of the value contained in
     *  the first element of periapsisTargetingInformation_ is given as output, because this element indicates whether
     *  the periapsis altitude falls within the boundaries of the corridor, whereas this function outputs whether the maneuver
     *  needs to be performed (which is by definition, the opposite).
     *  \return Boolean denoting whether the apoapsis maneuver is to be performed.
     */
    bool getIsApoapsisManeuverToBePerformed( ) { return !std::get< 0 >( periapsisTargetingInformation_ ); }

    //! Function to retrive the periapsis altitude targeting information.
    std::pair< double, double > getPeriapsisAltitudeTargetingInformation( )
    {
        return std::make_pair( std::get< 1 >( periapsisTargetingInformation_ ), std::get< 2 >( periapsisTargetingInformation_ ) );
    }

    //! Function to retirieve the value of the apoapsis maneuver vector.
    Eigen::Vector3d getScheduledApsisManeuver( ) { return scheduledApsisManeuver_; }

    //! Function to retrieve whether the input aerobraking phase is active.
    /*!
     *  Function to retrieve whether the the input aerobraking phase is active. The value of currentOrbitAerobrakingPhase_ is
     *  determined by the function determineAerobrakingPhase.
     *  \return Boolean denoting whether the input aerobraking phase is active.
     */
    bool getIsAerobrakingPhaseActive( const AerobrakingPhaseIndicator inputAerobrakingPhase )
    {
        return ( currentOrbitAerobrakingPhase_ == inputAerobrakingPhase );
    }

    //! Function to retrieve the history of estimated periapsis corridor boundaries.
    std::map< unsigned int, std::pair< double, double > > getHistoryOfEstimatedPeriapsisCorridorBoundaries( )
    {
        return historyOfEstimatedPeriapsisCorridorBoundaries_;
    }

    //! Function to retrieve the history of apo- and periapsis maneuver magnitudes.
    std::pair< double, std::map< unsigned int, double > > getHistoryOfApsisManeuverMagnitudes( )
    {
        // Sum all maneuver contributions
        double cumulativeManeuverMagnitude = 0;
        for ( std::map< unsigned int, double >::const_iterator mapIterator = historyOfApsisManeuverMagnitudes_.begin( );
              mapIterator != historyOfApsisManeuverMagnitudes_.end( ); mapIterator++ )
        {
            cumulativeManeuverMagnitude += std::fabs( mapIterator->second );
        }

        // Give output
        return std::make_pair( cumulativeManeuverMagnitude, historyOfApsisManeuverMagnitudes_ );
    }

private:

    //! Function to run corridor estimator with nominal conditions.
    /*!
     *  Function to run corridor estimator with nominal conditions.
     *  \param currentEstimatedCartesianState Vector denoting the current estimated translational Cartesian elements.
     *  \param currentEstimatedKeplerianState Vector denoting the current estimated translational Keplerian elements.
     *  \param useHeatAsLowerBoundaryThreshold Boolean denoting whether the heating conditions should be used as threshold,
     *      otherwise lifetime is used.
     */
    void estimateCorridorBoundaries( const double currentTime,
                                     const Eigen::Vector6d& currentEstimatedCartesianState,
                                     const Eigen::Vector6d& currentEstimatedKeplerianState,
                                     const bool useHeatAsLowerBoundaryThreshold = true );

    //! Function to compute the rotation from local to intertial frame.
    /*!
     *  Function to compute the rotation from local to intertial frame. First, the direction cosine matrix (DCM) describing the
     *  rotation from trajectory to inertial frame is found, then the rotation from local to trajectory is added. The first DCM is
     *  found by using the velocity and radial distance vector. The velocity (unit) vector corresponds directly to the x-axis of the
     *  trajectory frame, whereas the z-axis is computed by subtracting from the radial distance (unit) vector its projection on the
     *  x-axis, and then inverting to find the unit vector that points toward the planet. Then, the y-axis is determined via the
     *  right-hand rule (i.e., with the cross product). Note that since the estimated state is used, the actual transformation can differ.
     *  \param currentEstimatedCartesianState Current estimated Cartesian state as provided by the navigation system.
     *  \return Direction cosine matrix representing the estimated rotation from local to inertial frame.
     */
    Eigen::Matrix3d computeCurrentRotationFromLocalToInertialFrame( const Eigen::Vector6d& currentEstimatedCartesianState )
    {
        // Declare direction cosine matrix
        Eigen::Matrix3d transformationFromTrajectoryToInertialFrame;

        // Find the trajectory x-axis unit vector
        Eigen::Vector3d xUnitVector = currentEstimatedCartesianState.segment( 3, 3 ).normalized( );
        transformationFromTrajectoryToInertialFrame.col( 0 ) = xUnitVector;

        // Find trajectory z-axis unit vector
        Eigen::Vector3d zUnitVector = currentEstimatedCartesianState.segment( 0, 3 ).normalized( );
        zUnitVector -= zUnitVector.dot( xUnitVector ) * xUnitVector;
        zUnitVector *= -1.0; // axis points toward planet
        transformationFromTrajectoryToInertialFrame.col( 2 ) = zUnitVector;

        // Find body-fixed y-axis unit vector
        transformationFromTrajectoryToInertialFrame.col( 1 ) = zUnitVector.cross( xUnitVector );

        // Compute rotation from local to trajectory frame
        Eigen::Matrix3d transformationFromLocalToTrajectoryFrame = Eigen::Matrix3d::Zero( );
        transformationFromLocalToTrajectoryFrame( 0, 1 ) = 1.0;
        transformationFromLocalToTrajectoryFrame( 1, 2 ) = -1.0;
        transformationFromLocalToTrajectoryFrame( 2, 0 ) = -1.0;

        // Give output by combining the two rotations
        return transformationFromTrajectoryToInertialFrame * transformationFromLocalToTrajectoryFrame;
    }

    //! Function to compute the current estimated apoapsis radius.
    /*!
     *  Function to compute the current estimated apoapsis radius, where the assumption of a Kepler orbit (i.e., an unperturbed
     *  orbit) is taken. The actual result is then going to differ from this very simple first order estimation.
     *  \param currentEstimatedKeplerianState Current estimated translational Keplerian elements.
     *  \return First order approximation of the estimated apoapsis radius.
     */
    double computeCurrentFirstOrderEstimatedApoapsisRadius( const Eigen::Vector6d& currentEstimatedKeplerianState )
    {
        return currentEstimatedKeplerianState[ 0 ] * ( 1.0 + currentEstimatedKeplerianState[ 1 ] );
    }

    //! Function to compute the current estimated periapsis radius.
    /*!
     *  Function to compute the current estimated periapsis radius, where the assumption of a Kepler orbit (i.e., an unperturbed
     *  orbit) is taken. The actual result is then going to differ from this very simple first order estimation.
     *  \param currentEstimatedKeplerianState Current estimated translational Keplerian elements.
     *  \return First order approximation of the estimated periapsis radius.
     */
    double computeCurrentFirstOrderEstimatedPeriapsisRadius( const Eigen::Vector6d& currentEstimatedKeplerianState )
    {
        return currentEstimatedKeplerianState[ 0 ] * ( 1.0 - currentEstimatedKeplerianState[ 1 ] );
    }

    //! Double denoting the target periapsis altitude to achieve at the end of aerobraking.
    const double targetPeriapsisAltitude_;

    //! Double denoting the target apoapsis altitude to achieve at the end of aerobraking.
    const double targetApoapsisAltitude_;

    //! Double denoting the maximum allowed heat flux that the spacecraft can endure.
    const double maximumAllowedHeatRate_;

    //! Double denoting the maximum allowed heat load that the spacecraft can endure.
    const double maximumAllowedHeatLoad_;

    //! Double denoting the miminum allowed dynamic pressure that the spacecraft should encounter.
    const double minimumAllowedDynamicPressure_;

    //! Double denoting the minimum allowed predicted lifetime in days.
    const double minimumAllowedLifetime_;

    //! Boolean denoting whether the guidance system is being tested.
    const bool guidanceTesting_;

    //! Function returning the multiplication factor to be used to correct for the difference between the Kepler orbit assumption
    //! of the corridor estimator.
    boost::function< double( const double ) > altitudeCorrectionFunction_;

    //! Pointer to corridor estimator object.
    boost::shared_ptr< CorridorEstimator > corridorEstimator_;

    //! Pointer to root-finder used to esimate the value of the apoapsis maneuver.
    boost::shared_ptr< root_finders::BisectionCore< double > > maneuverBisectionRootFinder_;

    //! Standard gravitational parameter of body being orbited.
    double planetaryGravitationalParameter_;

    //! Radius of body being orbited.
    double planetaryRadius_;

    //! Function to propagate state with custom propagation termination settings.
    /*!
     *  Function to propagate state with custom propagation termination settings. This function simply refers to the
     *  propagateTranslationalStateWithCustomTerminationSettings function of the navigation system, in order to propagate the state
     *  from the input initial condition, until the custom termination conditions are reached.
     */
    boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > >(
            const boost::shared_ptr< propagators::PropagationTerminationSettings >, const Eigen::Vector6d& ) > statePropagationFunction_;

    //! Function to propagate state for two thirds of the orbital period.
    boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > >(
            const Eigen::Vector6d& ) > periodReducedStatePropagationFunction_;

    //! Function to propagate state for two times the lifetime limit.
    boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > >(
            const Eigen::Vector6d& ) > lifetimeReducedStatePropagationFunction_;

    //! Enumeration denoting the aerobraking phase belonging to the current orbit.
    AerobrakingPhaseIndicator currentOrbitAerobrakingPhase_;

    //! Double denoting the scaling to be applied to the periapsis altitude target altitude during the walk-in and walk-out phases.
    double periapsisAltitudeScaling_;

    //! Integer denoting the value of the current orbit.
    unsigned int currentOrbitCounter_;

    //! Tuple containing the periapsis targeting information.
    /*!
     *  Tuple containing the periapsis targeting information, where the first element is a boolean denoting whether the predicted
     *  periapsis altitude falls within the boundaries of the corridor, the second one denotes the predicted periapsis altitude,
     *  whereas the third and last element is the target periapsis altitude. If the first element is true, then the periapsis
     *  altitude is nothing but the predicted periapsis altitude.
     */
    std::tuple< bool, double, double > periapsisTargetingInformation_;

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsisManeuver_;

    //! History of estimated periapsis corridor boundaries.
    /*!
     *  History of estimated periapsis corridor boundaries, stored as a map, where the keys are the orbit numbers and the
     *  mapped values are pairs, whose first element is the lower altitude bound and the second one the upper altitude bound.
     */
    std::map< unsigned int, std::pair< double, double > > historyOfEstimatedPeriapsisCorridorBoundaries_;

    //! History of estimated apo- and periapsis maneuver magnitudes.
    std::map< unsigned int, double > historyOfApsisManeuverMagnitudes_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_GUIDANCE_SYSTEM_H
