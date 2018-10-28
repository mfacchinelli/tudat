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

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/SystemModels/UnitTests/generators.h"

using namespace tudat::system_models;

namespace tudat
{

namespace unit_tests
{

// Set test conditions
const std::pair< unsigned int, unsigned int > testConditions = { 0, 3 };
const std::pair< unsigned int, unsigned int > testModes = { 0, 2 };

BOOST_AUTO_TEST_SUITE( test_guidance_system )

//BOOST_AUTO_TEST_CASE( testGuidancePhase )
//{

//}

BOOST_AUTO_TEST_CASE( testCorridorEstimator )
{
    // Declare guidance system generator
    boost::shared_ptr< GuidanceSystemGenerator > systemGenerator = boost::make_shared< GuidanceSystemGenerator >( );
    double initialTime = systemGenerator->getCurrentTime( );
    double planetaryGravitationalParameter = systemGenerator->getGravitationalParameter( );
    double planetaryRadius = systemGenerator->getRadius( );

    // Declare guidance parameters
    const double targetPeriapsisAltitude = 255.0e3;
    const double targetApoapsisAltitude = 320.0e3;
    const double maximumAllowedHeatRate = 2800.0;
    const double maximumAllowedHeatLoad = 500.0e3;
    const double minimumAllowedDynamicPressure = 0.19;
    const double minimumAllowedLifetime = 2.0;

    // Generate guidance system
    boost::shared_ptr< GuidanceSystem > guidanceSystem = systemGenerator->createGuidanceSystem(
                targetPeriapsisAltitude, targetApoapsisAltitude, maximumAllowedHeatRate, maximumAllowedHeatLoad,
                minimumAllowedDynamicPressure, minimumAllowedLifetime );

    // Loop over various initial conditions
    for ( unsigned int initialConditions = testConditions.first; initialConditions < testConditions.second; initialConditions++ )
    {
        // Set initial conditions
        Eigen::Vector6d initialKeplerianState;
        std::pair< unsigned int, unsigned > pairOfAtmosphereInitiationIndicators;
        switch ( initialConditions )
        {
        case 0:
        {
            initialKeplerianState( 0 ) = 25946932.3;
            initialKeplerianState( 1 ) = 0.8651912;
            pairOfAtmosphereInitiationIndicators = { 1.0, 7.0 };
            break;
        }
        case 1:
        {
            initialKeplerianState( 0 ) = 6808709.3;
            initialKeplerianState( 1 ) = 0.4861342;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        case 2:
        {
            initialKeplerianState( 0 ) = 4699198.5;
            initialKeplerianState( 1 ) = 0.2546816;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        }
        initialKeplerianState( 2 ) = unit_conversions::convertDegreesToRadians( 93.0 );
        initialKeplerianState( 3 ) = unit_conversions::convertDegreesToRadians( 158.7 );
        initialKeplerianState( 4 ) = unit_conversions::convertDegreesToRadians( 23.4 );
        initialKeplerianState( 5 ) = unit_conversions::convertDegreesToRadians( 180.0 );

        // Convert to Cartesian elements
        Eigen::Vector6d initialCartesianState =
                orbital_element_conversions::convertKeplerianToCartesianElements( initialKeplerianState, planetaryGravitationalParameter );

        // Predict periapsis altitude with initial conditions
        unsigned int i = 0;
        double finalTime = initialTime + 2.0 / 3.0 *
                basic_astrodynamics::computeKeplerOrbitalPeriod( initialKeplerianState[ 0 ], planetaryGravitationalParameter );
        std::map< double, Eigen::VectorXd > trajectory = systemGenerator->propagateSpacecraftState( initialTime, finalTime,
                                                                                                    initialCartesianState ).second.first;
        Eigen::VectorXd historyOfAltitudes;
        historyOfAltitudes.resize( trajectory.size( ) );
        for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = trajectory.begin( );
              mapIterator != trajectory.end( ); mapIterator++, i++ )
        {
            historyOfAltitudes[ i ] = mapIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
        }
        double predictedPeriapsisAltitude = historyOfAltitudes.minCoeff( );

        // Loop over test modes
        for ( unsigned int navigationMode = testModes.first; navigationMode < testModes.second; navigationMode++ )
        {
            // Inform user
            std::cout << std::endl << "========================================================" << std::endl
                      << "Initial conditions: " << initialConditions << ". Test case: " << navigationMode << "." << std::endl
                      << "========================================================" << std::endl
                      << "Predicted periapsis altitude: " << predictedPeriapsisAltitude << std::endl;

            // Set estimated initial conditions
            Eigen::Vector6d estimatedCartesianState;
            switch ( navigationMode )
            {
            case 0:
            {
                // Use exact conditions
                estimatedCartesianState = initialCartesianState;
                break;
            }
            case 1:
            {
                // Generate noise distributions
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 50.0 }, 1 );
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 0.1 }, 2 );

                // Add noise to exact conditions
                estimatedCartesianState = initialCartesianState;
                for ( unsigned int i = 0; i < 6; i++ )
                {
                    if ( i < 3 )
                    {
                        estimatedCartesianState[ i ] += positionNoiseGenerator->getRandomVariableValue( );
                    }
                    else
                    {
                        estimatedCartesianState[ i ] += velocityNoiseGenerator->getRandomVariableValue( );
                    }
                }
                break;
            }
            }

            // Convert to Cartesian elements
            Eigen::Vector6d estimatedKeplerianState =
                    orbital_element_conversions::convertCartesianToKeplerianElements( estimatedCartesianState, planetaryGravitationalParameter );

            // Determine and check aerobraking phase
            guidanceSystem->determineAerobrakingPhase( estimatedKeplerianState, pairOfAtmosphereInitiationIndicators );
            switch ( initialConditions )
            {
            case 0:
                BOOST_CHECK_BITWISE_EQUAL( true, guidanceSystem->getIsAerobrakingPhaseActive( GuidanceSystem::walk_in_phase ) );
                break;
            case 1:
                BOOST_CHECK_BITWISE_EQUAL( true, guidanceSystem->getIsAerobrakingPhaseActive( GuidanceSystem::main_phase ) );
                break;
            case 2:
                BOOST_CHECK_BITWISE_EQUAL( true, guidanceSystem->getIsAerobrakingPhaseActive( GuidanceSystem::walk_out_phase ) );
                break;
            }

            // Run corridor estimator
            guidanceSystem->runCorridorEstimator( initialTime, estimatedCartesianState, estimatedKeplerianState );

            // Check estimated periapsis altitude
            BOOST_CHECK_CLOSE_FRACTION( predictedPeriapsisAltitude, guidanceSystem->getPeriapsisAltitudeTargetingInformation( ).first, 1.0e-2 );

            // Check that heating does not exceed maximum allowed values
            Eigen::Vector6d postManeuverCartesianState = initialCartesianState;
            guidanceSystem->runApoapsisManeuverEstimator(
                        estimatedCartesianState, estimatedKeplerianState,
                        basic_astrodynamics::computeKeplerMeanMotion( estimatedKeplerianState[ 0 ], planetaryGravitationalParameter ), true );
            postManeuverCartesianState.segment( 3, 3 ) += guidanceSystem->getScheduledApsisManeuver( );

            std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationData =
                    systemGenerator->propagateSpacecraftState( initialTime, finalTime, postManeuverCartesianState ).second;

            i = 0;
            Eigen::VectorXd historyOfTimes;
            Eigen::VectorXd historyOfAltitudes;
            Eigen::VectorXd historyOfDynamicPressure;
            Eigen::VectorXd historyOfHeatRate;
            historyOfTimes.resize( propagationData.first.size( ) );
            historyOfAltitudes.resize( propagationData.first.size( ) );
            historyOfDynamicPressure.resize( propagationData.first.size( ) );
            historyOfHeatRate.resize( propagationData.first.size( ) );
            for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = propagationData.first.begin( );
                  mapIterator != propagationData.first.end( ); mapIterator++, i++ )
            {
                historyOfTimes[ i ] = mapIterator->first;
                historyOfAltitudes[ i ] = mapIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
                historyOfDynamicPressure[ i ] = propagationData.second[ mapIterator->first ][ 0 ];
                historyOfHeatRate[ i ] = propagationData.second[ mapIterator->first ][ 1 ];
            }
            double postManeuverPredictedPeriapsisAltitude = historyOfAltitudes.minCoeff( );
            double peakDynamicPressure = historyOfDynamicPressure.maxCoeff( );
            double peakHeatRate = historyOfHeatRate.maxCoeff( );
            double heatLoad = numerical_quadrature::performTrapezoidalQuadrature(
                        utilities::convertEigenVectorToStlVector( historyOfTimes ),
                        utilities::convertEigenVectorToStlVector( historyOfHeatRate ) );
            std::cout << "Post-maneuver predicted periapsis altitude: " << postManeuverPredictedPeriapsisAltitude << std::endl
                      << "Peak dynamic pressure: " << peakDynamicPressure << " N/m^2" << std::endl
                      << "Peak heat rate: " << peakHeatRate << " W/m^2" << std::endl
                      << "Heat load: " << heatLoad / 1.0e3 << "kJ/m^2" << std::endl;

            // Check that target altitude was indeed targeted
            BOOST_CHECK_CLOSE_FRACTION( postManeuverPredictedPeriapsisAltitude,
                                        guidanceSystem->getPeriapsisAltitudeTargetingInformation( ).second, 1.0e-2 );

            // Check that heating conditions are met
            if ( !guidanceSystem->getIsAerobrakingPhaseActive( GuidanceSystem::walk_in_phase ) )
            {
                BOOST_CHECK_LE( minimumAllowedDynamicPressure, peakDynamicPressure );
            }
            BOOST_CHECK_GE( maximumAllowedHeatRate, peakHeatRate );
            BOOST_CHECK_GE( maximumAllowedHeatLoad, heatLoad );
        }
    }
}

//BOOST_AUTO_TEST_CASE( testManeuverEstimator )
//{

//}

BOOST_AUTO_TEST_SUITE_END( )

}

}

