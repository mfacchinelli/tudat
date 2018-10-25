/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

#include "Tudat/Astrodynamics/SystemModels/controlSystem.h"

using namespace tudat::system_models;

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_control_system )

BOOST_AUTO_TEST_CASE( testAttitudeController )
{
    // Define gains
    Eigen::Vector3d proportionalGain = Eigen::Vector3d::Constant( 0.75 );
    Eigen::Vector3d integralGain = Eigen::Vector3d::Constant( 0.3 );
    Eigen::Vector3d derivativeGain = Eigen::Vector3d::Constant( 1.5e-2 );

    // Declare control system
    boost::shared_ptr< ControlSystem > controlSystem = boost::make_shared< ControlSystem >( proportionalGain, integralGain, derivativeGain );

    // Check that attitude control vector is zero (since attitude control not implemented in this case)
    BOOST_CHECK_EQUAL( Eigen::Vector3d::Zero( ), controlSystem->getCurrentAttitudeControlVector( ) );
}

BOOST_AUTO_TEST_CASE( testOrbitController )
{
    // Define gains
    Eigen::Vector3d proportionalGain = Eigen::Vector3d::Constant( 0.75 );
    Eigen::Vector3d integralGain = Eigen::Vector3d::Constant( 0.3 );
    Eigen::Vector3d derivativeGain = Eigen::Vector3d::Constant( 1.5e-2 );

    // Declare control system
    boost::shared_ptr< ControlSystem > controlSystem = boost::make_shared< ControlSystem >( proportionalGain, integralGain, derivativeGain );

    // Check that if called, maneuver is zero (since it has not been set yet)
    BOOST_CHECK_EQUAL( Eigen::Vector3d::Zero( ), controlSystem->getScheduledApsisManeuver( ) );

    // Set apsis maneuver vector
    Eigen::Vector3d apsisManeuver = Eigen::Vector3d::Random( );
    controlSystem->updateOrbitController( apsisManeuver );

    // Check that maneuver is retrieved correctly
    BOOST_CHECK_EQUAL( apsisManeuver, controlSystem->getScheduledApsisManeuver( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

