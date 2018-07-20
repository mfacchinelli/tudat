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
#include <vector>

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Mathematics/RootFinders/bisection.h"

namespace tudat
{

namespace guidance_navigation_control
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
        aerobraking_completed = 3
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
                    const double minimumAllowedDynamicPressure ) :
        targetPeriapsisAltitude_( targetPeriapsisAltitude ), targetApoapsisAltitude_( targetApoapsisAltitude ),
        maximumAllowedHeatRate_( maximumAllowedHeatRate ), maximumAllowedHeatLoad_( maximumAllowedHeatLoad ),
        minimumAllowedDynamicPressure_( minimumAllowedDynamicPressure )
    {
        // Create root-finder object for bisection of periapsis altitude
        // The values inserted are the tolerance in independent value (i.e., the percentage corresponding to 0.5 km difference at
        // 100 km altitude) and the maximum number of iterations (i.e., 10 iterations)
        altitudeBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 0.5 / 100.0, 10 );
    }

    //! Destructor.
    ~GuidanceSystem( ) { }

    //! Function to determine in which aerobraking phase the spacecraft is currently in.
    /*!
     *  Function to determine in which aerobraking phase the spacecraft is currently in.
     *  \param currentEstimatedKeplerianState
     *  \param pairOfAtmosphereInitiationIndicators
     */
    void determineAerobrakingPhase( const Eigen::Vector6d& currentEstimatedKeplerianState,
                                    const std::pair< unsigned int, unsigned int >& pairOfAtmosphereInitiationIndicators )
    {
        // Declare aerobraking phase indicator and set value to main phase
        AerobrakingPhaseIndicator detectedAerobrakingPhase = main_phase;

        // Check whether atmosphere has been initiated
        if ( pairOfAtmosphereInitiationIndicators.first < pairOfAtmosphereInitiationIndicators.second )
        {
            detectedAerobrakingPhase = walk_in_phase;
            periapsisAltitudeWalkInScaling_ = 1.0;
        }

        // Check whether apoapsis is approaching target value
        double predictedApoapsisRadius = computeCurrentEstimatedApoapsisRadius( currentEstimatedKeplerianState );
        if ( std::fabs( predictedApoapsisRadius - targetApoapsisAltitude_ ) < 50e3 )
        {
            detectedAerobrakingPhase = walk_out_phase;
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
     *  \param currentEstimatedKeplerianState
     *  \param planetaryRadius
     *  \param planetaryGravitationalParameter
     *  \param atmosphericInterfaceRadius
     *  \param densityFunction
     */
    void runCorridorEstimator( const Eigen::Vector6d& currentEstimatedKeplerianState,
                               const double planetaryRadius,
                               const double planetaryGravitationalParameter,
                               const double atmosphericInterfaceRadius,
                               const boost::function< double( double ) >& densityFunction );

    //! Function to run maneuver estimator (ME).
    void runManeuverEstimator( const Eigen::Vector6d& currentEstimatedKeplerianState,
                               const double currentEstimatedMeanAnomaly,
                               const double planetaryRadius )
    {
        // Inform user
        std::cout << "Estimating Apoapsis Maneuver." << std::endl;

        // Set apoapsis maneuver vector to zero
        scheduledApsoapsisManeuver_.setZero( );

        // Compute predicted periapsis radius
        double predictedPeriapsisAltitude = computeCurrentEstimatedPeriapsisRadius( currentEstimatedKeplerianState ) - planetaryRadius;
        double differenceInPeriapsisAltitude = pairOfPeriapsisTargetingInformation_.second - predictedPeriapsisAltitude;

        // Compute estimated maneuver in y-direction of local orbit frame
        scheduledApsoapsisManeuver_[ 1 ] = 0.25 * currentEstimatedMeanAnomaly * differenceInPeriapsisAltitude * std::sqrt(
                    ( 1.0 + currentEstimatedKeplerianState[ 1 ] ) / ( 1.0 - currentEstimatedKeplerianState[ 1 ] ) );
        std::cout << "Scheduled maneuver: " << scheduledApsoapsisManeuver_.transpose( ) << std::endl;
    }

    //! Function to retrieve whether the apoapsis maneuver is to be performed.
    /*!
     *  Function to retrieve whether the apoapsis maneuver is to be performed. Note that the inverse of the value contained in
     *  the first element of pairOfPeriapsisTargetingInformation_ is given as output, because this element indicates whether
     *  the periapsis altitude falls within the boundaries of the corridor, whereas this function outputs whether the maneuver
     *  needs to be performed (which is by definition, the opposite).
     *  \return Boolean denoting whether the apoapsis maneuver is to be performed.
     */
    bool getIsApoapsisManeuverToBePerformed( ) { return !pairOfPeriapsisTargetingInformation_.first; }

    //! Function to retirieve the value of the apoapsis maneuver vector.
    Eigen::Vector3d getScheduledApoapsisManeuver( ) { return scheduledApsoapsisManeuver_; }

private:

    //! Function to compute the current estimated apoapsis radius.
    double computeCurrentEstimatedApoapsisRadius( const Eigen::Vector6d& currentEstimatedKeplerianState )
    {
        return currentEstimatedKeplerianState[ 0 ] * ( 1.0 + currentEstimatedKeplerianState[ 1 ] );
    }

    //! Function to compute the current estimated periapsis radius.
    double computeCurrentEstimatedPeriapsisRadius( const Eigen::Vector6d& currentEstimatedKeplerianState )
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

    //! Pointer to root-finder used to esimated the periapsis corridor altitudes.
    boost::shared_ptr< root_finders::BisectionCore< double > > altitudeBisectionRootFinder_;

    //! Enumeration denoting the aerobraking phase belonging to the current orbit.
    AerobrakingPhaseIndicator currentOrbitAerobrakingPhase_;

    //! Double denoting the scaling to be applied to the periapsis altitude target altitude during walk-in phase.
    double periapsisAltitudeWalkInScaling_;

    //! History of estimated periapsis corridor boundaries.
    /*!
     *  History of estimated periapsis corridor boundaries, stored as a map, where the keys are the orbit numbers and the
     *  mapped values are pairs, whose first element is the lower altitude bound and the second one the upper altitude bound.
     */
    std::vector< std::pair< double, double > > historyOfEstimatedPeriapsisCorridorBoundaries_;

    //! Pair containing the periapsis targeting information.
    /*!
     *  Pair containing the periapsis targeting information, where the first element is a boolean denoting whether the predicted
     *  periapsis altitude falls within the boundaries of the corridor, and the second one denotes the target periapsis altitude.
     *  If the first element is true, then the periapsis altitude is nothing but the predicted periapsis altitude.
     */
    std::pair< bool, double > pairOfPeriapsisTargetingInformation_;

    //! Vector denoting the velocity change scheduled to be applied at apoapsis.
    Eigen::Vector3d scheduledApsoapsisManeuver_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_GUIDANCE_SYSTEM_H
