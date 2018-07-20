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

    //! Constructor.
    GuidanceSystem( const double maximumAllowedHeatFlux,
                    const double maximumAllowedHeatLoad,
                    const double minimumAllowedDynamicPressure ) :
    maximumAllowedHeatFlux_( maximumAllowedHeatFlux ), maximumAllowedHeatLoad_( maximumAllowedHeatLoad ),
    minimumAllowedDynamicPressure_( minimumAllowedDynamicPressure )
    {
        // Create root-finder object for bisection of periapsis altitude
        // The values inserted are the tolerance in independent value (i.e., 0.5 km altitude) and the maximum
        // number of iterations (i.e., 10 iterations)
        altitudeBisectionRootFinder_ = boost::make_shared< root_finders::BisectionCore< double > >( 500, 10 );
    }

    //! Destructor.
    ~GuidanceSystem( ) { }

    //! Function to run corridor estimator (CE).
    void runCorridorEstimator( const unsigned int currentOrbitCounter );

    //! Function to run maneuver estimator (ME).
    void runManeuverEstimator( )
    {
        scheduledApsoapsisManeuver_ = Eigen::Vector3d::Zero( );
    }

    //! Function to retirieve the apoapsis maneuver.
    Eigen::Vector3d getScheduledApoapsisManeuver( ) { return scheduledApsoapsisManeuver_; }

private:

    //! Double denoting the maximum allowed heat flux that the spacecraft can endure.
    const double maximumAllowedHeatFlux_;

    //! Double denoting the maximum allowed heat load that the spacecraft can endure.
    const double maximumAllowedHeatLoad_;

    //! Double denoting the miminum allowed dynamic pressure that the spacecraft should encounter.
    const double minimumAllowedDynamicPressure_;

    //! Pointer to root-finder used to esimated the periapsis corridor altitudes.
    boost::shared_ptr< root_finders::BisectionCore< double > > altitudeBisectionRootFinder_;

    //! History of estimated periapsis corridor boundaries.
    /*!
     *  History of estimated periapsis corridor boundaries, stored as a map, where the keys are the orbit numbers and the
     *  mapped values are pairs, whose first element is the lower altitude bound and the second one the upper altitude bound.
     */
    std::map< int, std::pair< double, double > > historyOfEstimatedPeriapsisCorridorBoundaries_;

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
