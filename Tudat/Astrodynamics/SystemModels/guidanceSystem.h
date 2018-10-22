/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GUIDANCE_SYSTEM_H
#define TUDAT_GUIDANCE_SYSTEM_H

#include <map>
#include <tuple>
#include <vector>

#include <Eigen/Geometry>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Basics/utilities.h"

#include "Tudat/Mathematics/RootFinders/bisection.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"

#include "Tudat/Astrodynamics/Propagators/imanRmsMethod.h"

namespace tudat
{

namespace system_models
{

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
        aerobraking_complete = 3
    };

    //! Constructor.
    /*!
     *  Constructor.
     *  \param targetPeriapsisAltitude
     *  \param targetApoapsisAltitude
     *  \param maximumAllowedHeatRate
     *  \param maximumAllowedHeatLoad
     *  \param minimumAllowedDynamicPressure
     */
    GuidanceSystem( const double targetPeriapsisAltitude,
                    const double targetApoapsisAltitude,
                    const double maximumAllowedHeatRate,
                    const double maximumAllowedHeatLoad,
                    const double minimumAllowedDynamicPressure,
                    const double minimumAllowedLifetime ) :
        targetPeriapsisAltitude_( targetPeriapsisAltitude ), targetApoapsisAltitude_( targetApoapsisAltitude ),
        maximumAllowedHeatRate_( maximumAllowedHeatRate ), maximumAllowedHeatLoad_( maximumAllowedHeatLoad ),
        minimumAllowedDynamicPressure_( minimumAllowedDynamicPressure ), minimumAllowedLifetime_( minimumAllowedLifetime )
    {
        // Create root-finder object for bisection of periapsis altitude
        // The values inserted are the tolerance in independent value (i.e., the percentage corresponding to 0.5 km difference at
        // 100 km altitude) and the maximum number of iterations (i.e., 10 iterations)
        altitudeBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 0.5 / 100.0, 10 );

        // Create root-finder object for bisection of maneuver magnitude estimate
        // The values inserted are the tolerance in independent value (i.e., the percentage corresponding to 0.5 km difference at
        // 100 km altitude) and the maximum number of iterations (i.e., 10 iterations)
        maneuverBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 0.5 / 100.0, 10 );
    }

    //! Destructor.
    ~GuidanceSystem( ) { }

    //! Function to create the guidance system objects.
    /*!
     *  Function to create the guidance system objects.
     *  \param statePropagationFunction
     */
    void createGuidanceSystemObjects(
            const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
            std::map< double, Eigen::VectorXd > > >( const boost::shared_ptr< propagators::PropagationTerminationSettings >,
                const Eigen::Vector6d& ) >& statePropagationFunction )
    {
        // Create propagation function
        statePropagationFunction_ = statePropagationFunction;
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
        // Declare aerobraking phase indicator and set value to main phase
        AerobrakingPhaseIndicator detectedAerobrakingPhase = main_phase;

        // Set periapsis altitude scaling to default value
        periapsisAltitudeScaling_ = 1.0;

        // Check whether atmosphere has been initiated
        if ( pairOfAtmosphereInitiationIndicators.first < pairOfAtmosphereInitiationIndicators.second )
        {
            detectedAerobrakingPhase = walk_in_phase;
            periapsisAltitudeScaling_ = 1.2 - 0.2 *
                    static_cast< double >( pairOfAtmosphereInitiationIndicators.first ) /
                    static_cast< double >( pairOfAtmosphereInitiationIndicators.second );
        }

        // Check whether apoapsis is approaching target value
        double predictedApoapsisRadius = computeCurrentFirstOrderEstimatedApoapsisRadius( currentEstimatedKeplerianState );
        if ( ( std::fabs( predictedApoapsisRadius - targetApoapsisAltitude_ ) < 1500.0e3 ) ||
             ( currentEstimatedKeplerianState[ 1 ] < 0.3 ) )
        {
            detectedAerobrakingPhase = walk_out_phase;
        }
        if ( propagators::IMAN_ANALYSIS_INDEX != 0 )
        {
            std::cerr << "Periapsis altitude scaling OFF." << std::endl;
            periapsisAltitudeScaling_ = 1.0;
        }

        // Check whether aerobraking is complete
        double predictedPeriapsisRadius = computeCurrentFirstOrderEstimatedPeriapsisRadius( currentEstimatedKeplerianState );
        if ( ( std::fabs( predictedApoapsisRadius - targetApoapsisAltitude_ ) < 5.0e3 ) &&
             ( std::fabs( predictedPeriapsisRadius - targetPeriapsisAltitude_ ) < 5.0e3 ) )
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
     *  \param planetaryRadius Double denoting the radius of the planet being orbited.
     *  \param planetaryGravitationalParameter Double denoting the gravitational parameter of the planet being orbited.
     */
    void runCorridorEstimator( const double currentTime,
                               const Eigen::Vector6d& currentEstimatedCartesianState,
                               const Eigen::Vector6d& currentEstimatedKeplerianState,
                               const double planetaryRadius,
                               const double planetaryGravitationalParameter );

    //! Function to run maneuver estimator (ME).
    /*!
     *  Function to run maneuver estimator (ME). Note that since the root-finder uses the propagator defined in the runCorridorEstimator
     *  function (i.e., periodReducedStatePropagationFunction_), said function needs to be run before this one.
     *  \param currentEstimatedCartesianState Vector denoting the current estimated translational Cartesian elements.
     *  \param currentEstimatedKeplerianState Vector denoting the current estimated translational Keplerian elements.
     *  \param currentEstimatedMeanMotion Double denoting the current estimated mean motion.
     *  \param planetaryRadius Double denoting the radius of the planet being orbited.
     *  \param improveEstimateWithBisection Boolean denoting whether the maneuver estimate should be improved by using a
     *      bisection root-finder algorithm.
     */
    void runManeuverEstimator( const Eigen::Vector6d& currentEstimatedCartesianState,
                               const Eigen::Vector6d& currentEstimatedKeplerianState,
                               const double currentEstimatedMeanMotion,
                               const double planetaryRadius,
                               const bool improveEstimateWithBisection = true );

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

    //! Function to retirieve the value of the apoapsis maneuver vector.
    Eigen::Vector3d getScheduledApoapsisManeuver( ) { return scheduledApsoapsisManeuver_; }

    //! Function to retrieve whether the aerobraking maneuver has been completed.
    /*!
     *  Function to retrieve whether the aerobraking maneuver has been completed. The value of currentOrbitAerobrakingPhase_ is
     *  determined by the function determineAerobrakingPhase and returns that the aerobraking phase has been completed, if the
     *  estimated periapsis and apoapsis altitudes are within 5 km of their respective target values.
     *  \return Boolean denoting whether the aerobraking maneuver has been completed.
     */
    bool getIsAerobrakingComplete( ) { return ( currentOrbitAerobrakingPhase_ == aerobraking_complete ); }

    //! Function to retrieve the history of estimated periapsis corridor boundaries.
    std::map< unsigned int, std::pair< double, double > > getHistoryOfEstimatedPeriapsisCorridorBoundaries( )
    {
        return historyOfEstimatedPeriapsisCorridorBoundaries_;
    }

    //! Function to retrieve the history of apoapsis maneuver magnitudes.
    std::pair< double, std::map< unsigned int, double > > getHistoryOfApoapsisManeuverMagnitudes( )
    {
        // Sum all maneuver contributions
        double cumulativeManeuverMagnitude = 0;
        for ( std::map< unsigned int, double >::const_iterator mapIterator = historyOfApoapsisManeuverMagnitudes_.begin( );
              mapIterator != historyOfApoapsisManeuverMagnitudes_.end( ); mapIterator++ )
        {
            cumulativeManeuverMagnitude += mapIterator->second;
        }

        // Give output
        return std::make_pair( cumulativeManeuverMagnitude, historyOfApoapsisManeuverMagnitudes_ );
    }

private:

    //! Function to run corridor estimator with nominal conditions.
    /*!
     *  Function to run corridor estimator with nominal conditions.
     *  \param currentEstimatedCartesianState Vector denoting the current estimated translational Cartesian elements.
     *  \param currentEstimatedKeplerianState Vector denoting the current estimated translational Keplerian elements.
     *  \param planetaryRadius Double denoting the radius of the planet being orbited.
     *  \param planetaryGravitationalParameter Double denoting the gravitational parameter of the planet being orbited.
     *  \param useHeatAsLowerBoundaryThreshold Boolean denoting whether the heating conditions should be used as threshold,
     *      otherwise lifetime is used.
     */
    void estimateCorridorBoundaries( const double currentTime,
                                     const Eigen::Vector6d& currentEstimatedCartesianState,
                                     const Eigen::Vector6d& currentEstimatedKeplerianState,
                                     const double planetaryRadius,
                                     const double planetaryGravitationalParameter,
                                     const bool useHeatAsLowerBoundaryThreshold = true );

    //! Function to compute the rotation from local to intertial frame.
    /*!
     *  Function to compute the rotation from local to intertial frame. First, the direction cosine matrix (DCM) describing the
     *  rotation from trajectory to inertial frame is found, then the rotation from local to trajectory is added. The first DCM is
     *  found by using the velocity and radial distance vector. The velocity (unit) vector corresponds directly to the x-axis of the
     *  trajectory frame, whereas the z-axis is computed by subtracting from the radial distance (unit) vector its projection on the
     *  x-axis. Then, the y-axis is determined via the right-hand rule (i.e., with the cross product). Note that since the estimated
     *  state is used, the actual transformation can differ.
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

    //! Pointer to root-finder used to esimate the periapsis corridor altitudes.
    boost::shared_ptr< root_finders::BisectionCore< double > > altitudeBisectionRootFinder_;

    //! Pointer to root-finder used to esimate the value of the apoapsis maneuver.
    boost::shared_ptr< root_finders::BisectionCore< double > > maneuverBisectionRootFinder_;

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
    Eigen::Vector3d scheduledApsoapsisManeuver_;

    //! History of estimated periapsis corridor boundaries.
    /*!
     *  History of estimated periapsis corridor boundaries, stored as a map, where the keys are the orbit numbers and the
     *  mapped values are pairs, whose first element is the lower altitude bound and the second one the upper altitude bound.
     */
    std::map< unsigned int, std::pair< double, double > > historyOfEstimatedPeriapsisCorridorBoundaries_;

    //! History of estimated apoapsis maneuver magnitudes.
    std::map< unsigned int, double > historyOfApoapsisManeuverMagnitudes_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_GUIDANCE_SYSTEM_H
