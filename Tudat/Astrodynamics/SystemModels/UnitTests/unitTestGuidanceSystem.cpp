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
#include "Tudat/InputOutput/multiDimensionalArrayWriter.h"

using namespace tudat::system_models;

namespace tudat
{

namespace unit_tests
{

// Set test conditions
const std::pair< unsigned int, unsigned int > testConditions = { 0, 3 };
const std::pair< unsigned int, unsigned int > testModes = { 0, 2 };
const std::pair< int, int > testMagnitudes = { -5, 4 };
const std::pair< int, int > testAltitudes = { -6, 3 };
const unsigned int numberOfMonteCarlos = 20;

BOOST_AUTO_TEST_SUITE( test_guidance_system )

//BOOST_AUTO_TEST_CASE( testAerobrakingPhase )
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
    const double maximumAllowedHeatRate = 2800.0 / 2.0;
    const double maximumAllowedHeatLoad = 500.0e3 / 2.0;
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
            initialKeplerianState( 0 ) = 5287404.0;
            initialKeplerianState( 1 ) = 0.3385648;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        case 2:
        {
            initialKeplerianState( 0 ) = 4710917.5;
            initialKeplerianState( 1 ) = 0.2515604;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        }
        initialKeplerianState( 2 ) = unit_conversions::convertDegreesToRadians( 93.0 );
        initialKeplerianState( 4 ) = unit_conversions::convertDegreesToRadians( 158.7 );
        initialKeplerianState( 3 ) = unit_conversions::convertDegreesToRadians( 43.6 );
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
                      << "Predicted periapsis altitude: " << predictedPeriapsisAltitude / 1.0e3 << std::endl;

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

            // Run maneuver estimator
            Eigen::Vector6d postManeuverCartesianState = initialCartesianState;
            guidanceSystem->runApoapsisManeuverEstimator( estimatedCartesianState, estimatedKeplerianState, true );
            postManeuverCartesianState.segment( 3, 3 ) += guidanceSystem->getScheduledApsisManeuver( );

            // Propagate state
            std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationData =
                    systemGenerator->propagateSpacecraftState( initialTime, finalTime, postManeuverCartesianState ).second;

            // Check that heating does not exceed maximum allowed values
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
            std::cout << "Post-maneuver predicted periapsis altitude: " << postManeuverPredictedPeriapsisAltitude / 1.0e3 << std::endl
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

BOOST_AUTO_TEST_CASE( testManeuverEstimator )
{

    // Declare guidance system generator
    boost::shared_ptr< GuidanceSystemGenerator > systemGenerator = boost::make_shared< GuidanceSystemGenerator >( );
    double initialTime = systemGenerator->getCurrentTime( );
    double planetaryGravitationalParameter = systemGenerator->getGravitationalParameter( );
    double planetaryRadius = systemGenerator->getRadius( );

    // Declare guidance parameters
    const double targetFinalPeriapsisAltitude = 255.0e3;
    const double targetFinalApoapsisAltitude = 320.0e3;
    const double maximumAllowedHeatRate = 2800.0 / 2.0;
    const double maximumAllowedHeatLoad = 500.0e3 / 2.0;
    const double minimumAllowedDynamicPressure = 0.19;
    const double minimumAllowedLifetime = 2.0;

    // Generate guidance system
    boost::shared_ptr< GuidanceSystem > guidanceSystem = systemGenerator->createGuidanceSystem(
                targetFinalPeriapsisAltitude, targetFinalApoapsisAltitude, maximumAllowedHeatRate, maximumAllowedHeatLoad,
                minimumAllowedDynamicPressure, minimumAllowedLifetime );

    // Create multi-array of results
    const int magnitudeSpan = testMagnitudes.second - testMagnitudes.first;
    boost::multi_array< Eigen::Vector6d, 3 > aerothermodynamicsConditions(
                boost::extents[ testConditions.second ][ testModes.second ][ magnitudeSpan ] );

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
            initialKeplerianState( 0 ) = 5287404.0;
            initialKeplerianState( 1 ) = 0.3385648;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        case 2:
        {
            initialKeplerianState( 0 ) = 4710917.5;
            initialKeplerianState( 1 ) = 0.2515604;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        }
        initialKeplerianState( 2 ) = unit_conversions::convertDegreesToRadians( 93.0 );
        initialKeplerianState( 4 ) = unit_conversions::convertDegreesToRadians( 158.7 );
        initialKeplerianState( 3 ) = unit_conversions::convertDegreesToRadians( 43.6 );
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
                      << "Predicted periapsis altitude: " << predictedPeriapsisAltitude / 1.0e3 << std::endl;

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

            // Run corridor estimator
            guidanceSystem->runCorridorEstimator( initialTime, estimatedCartesianState, estimatedKeplerianState );
            double targetPeriapsisAltitude = guidanceSystem->getPeriapsisAltitudeTargetingInformation( ).second;

            // Run maneuver estimator
            guidanceSystem->runApoapsisManeuverEstimator( estimatedCartesianState, estimatedKeplerianState, true );
            double nominalManeuverMagnitude = guidanceSystem->getHistoryOfApsisManeuverMagnitudes( ).second.rbegin( )->second;
            Eigen::Matrix3d currentRotation = guidanceSystem->testComputeCurrentRotationFromLocalToInertialFrame( estimatedCartesianState );

            // Loop over various altitude differences
            for ( int magnitudeDifference = testMagnitudes.first; magnitudeDifference < testMagnitudes.second; magnitudeDifference++ )
            {
                // Reset maenuver magnitude
                double adjustedManeuverMagnitude = nominalManeuverMagnitude + 3.0 * static_cast< double >( magnitudeDifference ) * 1.0e-1;
                Eigen::Vector3d adjustedManeuver = Eigen::Vector3d::Zero( );
                adjustedManeuver[ 1 ] = adjustedManeuverMagnitude;
                adjustedManeuver = currentRotation * adjustedManeuver;
                Eigen::Vector6d postManeuverCartesianState = initialCartesianState;
                postManeuverCartesianState.segment( 3, 3 ) += adjustedManeuver;

                // Inform user
                std::cout << std::endl << "====================== "
                          << "Magnitude difference: " << magnitudeDifference
                          << " ======================" << std::endl
                          << "Nominal magnitude: " << nominalManeuverMagnitude << std::endl
                          << "Adjusted magnitude: " << adjustedManeuverMagnitude << std::endl;

                // Propagate state
                std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationData =
                        systemGenerator->propagateSpacecraftState( initialTime, finalTime, postManeuverCartesianState ).second;

                // Check that heating does not exceed maximum allowed values
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
                std::cout << "Target periapsis altitude: " << targetPeriapsisAltitude / 1.0e3 << std::endl
                          << "Post-maneuver predicted periapsis altitude: " << postManeuverPredictedPeriapsisAltitude / 1.0e3 << std::endl
                          << "Peak dynamic pressure: " << peakDynamicPressure << " N/m^2" << std::endl
                          << "Peak heat rate: " << peakHeatRate << " W/m^2" << std::endl
                          << "Heat load: " << heatLoad / 1.0e3 << " kJ/m^2" << std::endl;

                // Store results
                aerothermodynamicsConditions[ initialConditions ][ navigationMode ][ magnitudeDifference - testMagnitudes.first ] =
                        ( Eigen::VectorXd( 6 ) << peakDynamicPressure, peakHeatRate, heatLoad,
                          nominalManeuverMagnitude, targetPeriapsisAltitude, postManeuverPredictedPeriapsisAltitude ).finished( );

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

    // Save results to file
    std::map< int, std::string > outFileNamesMap;
    outFileNamesMap[ 0 ] = "TestingResults/gsMePeakDynPress.dat";
    outFileNamesMap[ 1 ] = "TestingResults/gsMePeakHeatRate.dat";
    outFileNamesMap[ 2 ] = "TestingResults/gsMeHeatLoad.dat";
    outFileNamesMap[ 3 ] = "TestingResults/gsMeNominalManeuver.dat";
    outFileNamesMap[ 4 ] = "TestingResults/gsMeTargetPeriapsis.dat";
    outFileNamesMap[ 5 ] = "TestingResults/gsMeActualPeriapsis.dat";
    Eigen::VectorXd vectorOfConditions = Eigen::ArrayXd::LinSpaced( ( testConditions.second - testConditions.first ),
                                                                    testConditions.first, testConditions.second - 1 );
    Eigen::VectorXd vectorOfModes = Eigen::ArrayXd::LinSpaced( ( testModes.second - testModes.first ),
                                                               testModes.first, testModes.second - 1 );
    Eigen::VectorXd vectorOfMagnitudes = Eigen::ArrayXd::LinSpaced( ( testMagnitudes.second - testMagnitudes.first ),
                                                                   testMagnitudes.first, testMagnitudes.second - 1 );
    std::vector< std::vector< double > > independentVariables;
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfConditions ) );
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfModes ) );
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfMagnitudes ) );
    input_output::MultiArrayFileWriter< 3, 6 >::writeMultiArrayAndIndependentVariablesToFiles( outFileNamesMap,
                                                                                               independentVariables,
                                                                                               aerothermodynamicsConditions );
}

BOOST_AUTO_TEST_CASE( testPeriapsisSensitivity )
{
    // Declare guidance system generator
    boost::shared_ptr< GuidanceSystemGenerator > systemGenerator = boost::make_shared< GuidanceSystemGenerator >( );
    double initialTime = systemGenerator->getCurrentTime( );
    double planetaryGravitationalParameter = systemGenerator->getGravitationalParameter( );
    double planetaryRadius = systemGenerator->getRadius( );

    // Declare guidance parameters
    const double targetFinalPeriapsisAltitude = 255.0e3;
    const double targetFinalApoapsisAltitude = 320.0e3;
    const double maximumAllowedHeatRate = 2800.0 / 2.0;
    const double maximumAllowedHeatLoad = 500.0e3 / 2.0;
    const double minimumAllowedDynamicPressure = 0.19;
    const double minimumAllowedLifetime = 2.0;

    // Generate guidance system
    boost::shared_ptr< GuidanceSystem > guidanceSystem = systemGenerator->createGuidanceSystem(
                targetFinalPeriapsisAltitude, targetFinalApoapsisAltitude, maximumAllowedHeatRate, maximumAllowedHeatLoad,
                minimumAllowedDynamicPressure, minimumAllowedLifetime );

    // Create multi-array of results
    const int altitudeSpan = testAltitudes.second - testAltitudes.first;
    boost::multi_array< Eigen::Vector6d, 3 > aerothermodynamicsConditions(
                boost::extents[ testConditions.second ][ testModes.second ][ altitudeSpan ] );

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
            initialKeplerianState( 0 ) = 5287404.0;
            initialKeplerianState( 1 ) = 0.3385648;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        case 2:
        {
            initialKeplerianState( 0 ) = 4710917.5;
            initialKeplerianState( 1 ) = 0.2515604;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        }
        initialKeplerianState( 2 ) = unit_conversions::convertDegreesToRadians( 93.0 );
        initialKeplerianState( 4 ) = unit_conversions::convertDegreesToRadians( 158.7 );
        initialKeplerianState( 3 ) = unit_conversions::convertDegreesToRadians( 43.6 );
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
                      << "Predicted periapsis altitude: " << predictedPeriapsisAltitude / 1.0e3 << std::endl;

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

            // Run corridor estimator
            guidanceSystem->runCorridorEstimator( initialTime, estimatedCartesianState, estimatedKeplerianState );
            std::pair< double, double > nominalPeriapsisTargetingInformation = guidanceSystem->getPeriapsisAltitudeTargetingInformation( );

            // Loop over various altitude differences
            for ( int altitudeDifference = testAltitudes.first; altitudeDifference < testAltitudes.second; altitudeDifference++ )
            {
                // Reset periapsis targeting information
                std::pair< double, double > adjustedPeriapsisTargetingInformation = nominalPeriapsisTargetingInformation;
                adjustedPeriapsisTargetingInformation.second += 1.5 * static_cast< double >( altitudeDifference ) * 1.0e3;
                guidanceSystem->testSetPeriapsisAltitudeTargetingInformation( adjustedPeriapsisTargetingInformation );

                // Inform user
                std::cout << std::endl << "====================== "
                          << "Altitude difference: " << altitudeDifference
                          << " ======================" << std::endl
                          << "Nominal target periapsis altitude: " << nominalPeriapsisTargetingInformation.second / 1.0e3 << std::endl
                          << "Adjusted target periapsis altitude: " << adjustedPeriapsisTargetingInformation.second / 1.0e3 << std::endl;

                // Run maneuver estimator
                Eigen::Vector6d postManeuverCartesianState = initialCartesianState;
                guidanceSystem->runApoapsisManeuverEstimator( estimatedCartesianState, estimatedKeplerianState, true );
                postManeuverCartesianState.segment( 3, 3 ) += guidanceSystem->getScheduledApsisManeuver( );

                // Propagate state
                std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationData =
                        systemGenerator->propagateSpacecraftState( initialTime, finalTime, postManeuverCartesianState ).second;

                // Check that heating does not exceed maximum allowed values
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
                std::cout << "Post-maneuver predicted periapsis altitude: " << postManeuverPredictedPeriapsisAltitude / 1.0e3 << std::endl
                          << "Peak dynamic pressure: " << peakDynamicPressure << " N/m^2" << std::endl
                          << "Peak heat rate: " << peakHeatRate << " W/m^2" << std::endl
                          << "Heat load: " << heatLoad / 1.0e3 << " kJ/m^2" << std::endl;

                // Store results
                aerothermodynamicsConditions[ initialConditions ][ navigationMode ][ altitudeDifference - testAltitudes.first ] =
                        ( Eigen::VectorXd( 6 ) << peakDynamicPressure, peakHeatRate, heatLoad,
                          nominalPeriapsisTargetingInformation.second, adjustedPeriapsisTargetingInformation.second,
                          postManeuverPredictedPeriapsisAltitude ).finished( );

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

    // Save results to file
    std::map< int, std::string > outFileNamesMap;
    outFileNamesMap[ 0 ] = "TestingResults/gsPsPeakDynPress.dat";
    outFileNamesMap[ 1 ] = "TestingResults/gsPsPeakHeatRate.dat";
    outFileNamesMap[ 2 ] = "TestingResults/gsPsHeatLoad.dat";
    outFileNamesMap[ 3 ] = "TestingResults/gsPsNominalTarget.dat";
    outFileNamesMap[ 4 ] = "TestingResults/gsPsAdjustedTarget.dat";
    outFileNamesMap[ 5 ] = "TestingResults/gsPsActualPeriapsis.dat";
    Eigen::VectorXd vectorOfConditions = Eigen::ArrayXd::LinSpaced( ( testConditions.second - testConditions.first ),
                                                                    testConditions.first, testConditions.second - 1 );
    Eigen::VectorXd vectorOfModes = Eigen::ArrayXd::LinSpaced( ( testModes.second - testModes.first ),
                                                               testModes.first, testModes.second - 1 );
    Eigen::VectorXd vectorOfAltitudes = Eigen::ArrayXd::LinSpaced( ( testAltitudes.second - testAltitudes.first ),
                                                                   testAltitudes.first, testAltitudes.second - 1 );
    std::vector< std::vector< double > > independentVariables;
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfConditions ) );
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfModes ) );
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfAltitudes ) );
    input_output::MultiArrayFileWriter< 3, 6 >::writeMultiArrayAndIndependentVariablesToFiles( outFileNamesMap,
                                                                                               independentVariables,
                                                                                               aerothermodynamicsConditions );
}

BOOST_AUTO_TEST_CASE( testPeriapsisMonteCarlo )
{
    // Declare guidance system generator
    boost::shared_ptr< GuidanceSystemGenerator > systemGenerator = boost::make_shared< GuidanceSystemGenerator >( );
    double initialTime = systemGenerator->getCurrentTime( );
    double planetaryGravitationalParameter = systemGenerator->getGravitationalParameter( );
    double planetaryRadius = systemGenerator->getRadius( );

    // Declare guidance parameters
    const double targetFinalPeriapsisAltitude = 255.0e3;
    const double targetFinalApoapsisAltitude = 320.0e3;
    const double maximumAllowedHeatRate = 2800.0 / 2.0;
    const double maximumAllowedHeatLoad = 500.0e3 / 2.0;
    const double minimumAllowedDynamicPressure = 0.19;
    const double minimumAllowedLifetime = 2.0;

    // Generate guidance system
    boost::shared_ptr< GuidanceSystem > guidanceSystem = systemGenerator->createGuidanceSystem(
                targetFinalPeriapsisAltitude, targetFinalApoapsisAltitude, maximumAllowedHeatRate, maximumAllowedHeatLoad,
                minimumAllowedDynamicPressure, minimumAllowedLifetime );

    // Create multi-array of results
    boost::multi_array< Eigen::Matrix< double, 11, 1 >, 2 > aerothermodynamicsConditions(
                boost::extents[ numberOfMonteCarlos ][ testConditions.second ] );

    // Loop over various initial conditions
    for ( unsigned int initialConditions = testConditions.first; initialConditions < testConditions.second; initialConditions++ )
    {
        // Inform user
        std::cout << std::endl << "========================================================" << std::endl
                  << "Initial conditions: " << initialConditions << std::endl
                  << "========================================================" << std::endl;

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
            initialKeplerianState( 0 ) = 5287404.0;
            initialKeplerianState( 1 ) = 0.3385648;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        case 2:
        {
            initialKeplerianState( 0 ) = 4710917.5;
            initialKeplerianState( 1 ) = 0.2515604;
            pairOfAtmosphereInitiationIndicators = { 10.0, 7.0 };
            break;
        }
        }
        initialKeplerianState( 2 ) = unit_conversions::convertDegreesToRadians( 93.0 );
        initialKeplerianState( 4 ) = unit_conversions::convertDegreesToRadians( 158.7 );
        initialKeplerianState( 3 ) = unit_conversions::convertDegreesToRadians( 43.6 );
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

        // Generate noise distributions
        boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                statistics::createBoostContinuousRandomVariableGenerator(
                    statistics::normal_boost_distribution, { 0.0, 50.0 }, 1 );
        boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                statistics::createBoostContinuousRandomVariableGenerator(
                    statistics::normal_boost_distribution, { 0.0, 0.1 }, 2 );

        // Determine nominal condition
        guidanceSystem->determineAerobrakingPhase( initialKeplerianState, pairOfAtmosphereInitiationIndicators );
        guidanceSystem->runCorridorEstimator( initialTime, initialCartesianState, initialKeplerianState );
        std::pair< double, double > nominalPeriapsisTargetingInformation = guidanceSystem->getPeriapsisAltitudeTargetingInformation( );

        // Loop over various altitude differences
        for ( unsigned int monteCarloSimulation = 0; monteCarloSimulation < numberOfMonteCarlos; monteCarloSimulation++ )
        {
            // Perturb initial state
            Eigen::Vector6d estimatedCartesianState = initialCartesianState;
            Eigen::Vector6d errorCartesianState;
            for ( unsigned int i = 0; i < 6; i++ )
            {
                if ( i < 3 )
                {
                    errorCartesianState[ i ] = positionNoiseGenerator->getRandomVariableValue( );
                }
                else
                {
                    errorCartesianState[ i ] = velocityNoiseGenerator->getRandomVariableValue( );
                }
                estimatedCartesianState[ i ] += errorCartesianState[ i ];
            }

            // Convert to Cartesian elements
            Eigen::Vector6d estimatedKeplerianState =
                    orbital_element_conversions::convertCartesianToKeplerianElements( estimatedCartesianState, planetaryGravitationalParameter );

            // Determine aerobraking phase
            guidanceSystem->determineAerobrakingPhase( estimatedKeplerianState, pairOfAtmosphereInitiationIndicators );

            // Run corridor estimator
            guidanceSystem->runCorridorEstimator( initialTime, estimatedCartesianState, estimatedKeplerianState );
            std::pair< double, double > adjustedPeriapsisTargetingInformation = guidanceSystem->getPeriapsisAltitudeTargetingInformation( );

            // Run maneuver estimator
            Eigen::Vector6d postManeuverCartesianState = initialCartesianState;
            guidanceSystem->runApoapsisManeuverEstimator( estimatedCartesianState, estimatedKeplerianState, true );
            postManeuverCartesianState.segment( 3, 3 ) += guidanceSystem->getScheduledApsisManeuver( );

            // Propagate state
            std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationData =
                    systemGenerator->propagateSpacecraftState( initialTime, finalTime, postManeuverCartesianState ).second;

            // Check that heating does not exceed maximum allowed values
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
            std::cout << "Post-maneuver predicted periapsis altitude: " << postManeuverPredictedPeriapsisAltitude / 1.0e3 << std::endl
                      << "Peak dynamic pressure: " << peakDynamicPressure << " N/m^2" << std::endl
                      << "Peak heat rate: " << peakHeatRate << " W/m^2" << std::endl
                      << "Heat load: " << heatLoad / 1.0e3 << " kJ/m^2" << std::endl;

            // Store results
            aerothermodynamicsConditions[ monteCarloSimulation ][ initialConditions ] =
                    ( Eigen::VectorXd( 11 ) << peakDynamicPressure, peakHeatRate, heatLoad,
                      nominalPeriapsisTargetingInformation.second, adjustedPeriapsisTargetingInformation.second,
                      errorCartesianState ).finished( );

            // Check that heating conditions are met
            if ( !guidanceSystem->getIsAerobrakingPhaseActive( GuidanceSystem::walk_in_phase ) )
            {
                BOOST_CHECK_LE( minimumAllowedDynamicPressure, peakDynamicPressure );
            }
            BOOST_CHECK_GE( maximumAllowedHeatRate, peakHeatRate );
            BOOST_CHECK_GE( maximumAllowedHeatLoad, heatLoad );
        }
    }

    // Save results to file
    std::map< int, std::string > outFileNamesMap;
    outFileNamesMap[ 0 ] = "TestingResults/gsMcPeakDynPress.dat";
    outFileNamesMap[ 1 ] = "TestingResults/gsMcPeakHeatRate.dat";
    outFileNamesMap[ 2 ] = "TestingResults/gsMcHeatLoad.dat";
    outFileNamesMap[ 3 ] = "TestingResults/gsMcNominalTarget.dat";
    outFileNamesMap[ 4 ] = "TestingResults/gsMcAdjustedTarget.dat";
    outFileNamesMap[ 5 ] = "TestingResults/gsMcXError.dat";
    outFileNamesMap[ 6 ] = "TestingResults/gsMcYError.dat";
    outFileNamesMap[ 7 ] = "TestingResults/gsMcZError.dat";
    outFileNamesMap[ 8 ] = "TestingResults/gsMcVxError.dat";
    outFileNamesMap[ 9 ] = "TestingResults/gsMcVyError.dat";
    outFileNamesMap[ 10 ] = "TestingResults/gsMcVzError.dat";
    Eigen::VectorXd vectorOfConditions = Eigen::ArrayXd::LinSpaced( ( testConditions.second - testConditions.first ),
                                                                    testConditions.first, testConditions.second - 1 );
    Eigen::VectorXd vectorOfSimulations = Eigen::ArrayXd::LinSpaced( numberOfMonteCarlos, 0, numberOfMonteCarlos - 1 );
    std::vector< std::vector< double > > independentVariables;
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfSimulations ) );
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfConditions ) );
    input_output::MultiArrayFileWriter< 2, 11 >::writeMultiArrayAndIndependentVariablesToFiles( outFileNamesMap,
                                                                                                independentVariables,
                                                                                                aerothermodynamicsConditions );
}

//BOOST_AUTO_TEST_CASE( testPeriapsisSensitivity )
//{

//}

BOOST_AUTO_TEST_SUITE_END( )

}

}

