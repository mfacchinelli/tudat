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

#include <functional>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/SystemModels/UnitTests/generators.h"

#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

using namespace tudat::system_models;

namespace tudat
{

namespace unit_tests
{

// Set test conditions
const std::pair< unsigned int, unsigned int > testAtmospheres = { 0, 3 };
const std::pair< unsigned int, unsigned int > testConditions = { 0, 3 };
const std::pair< unsigned int, unsigned int > testModes = { 0, 1 };

// Linear offset function
double offsetFunctionWrapper( const double initialTime, const double time, const double acceleration )
{
    return ( 5.0e-2 + 1e1 * acceleration ) * ( time - initialTime );
}

BOOST_AUTO_TEST_SUITE( test_navigation_system )

//BOOST_AUTO_TEST_CASE( testNavigationPhase )
//{

//}

BOOST_AUTO_TEST_CASE( testStateEstimator )
{
    // Declare navigation system generator
    boost::shared_ptr< NavigationSystemGenerator > systemGenerator =
            boost::make_shared< NavigationSystemGenerator >( aerodynamics::exponential_atmosphere_model );

    // Define constants
    const double onboardTimeStep = 0.1;

    // Create multi-array of results
    boost::multi_array< Eigen::Vector4d, 3 > propagationResults(
                boost::extents[ 2 ][ testConditions.second ][ testModes.second ] );

    // Loop over Sun gravity usage
    for ( unsigned int sunGravityIndex = 0; sunGravityIndex < 2; sunGravityIndex++ )
    {
        // Loop over various initial conditions
        for ( unsigned int initialConditions = testConditions.first; initialConditions < testConditions.second; initialConditions++ )
        {
            // Propagate spacecraft state and retireve solution
            const std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > fullNumericalResults =
                    propagateStateForAerobrakingScenario( initialConditions, static_cast< bool >( sunGravityIndex ) );
            const std::map< double, Eigen::VectorXd > numericalIntegrationResult = fullNumericalResults.first;
            const std::map< double, Eigen::VectorXd > dependentVariablesResult = fullNumericalResults.second;

            // Extract times
            const double initialTime = numericalIntegrationResult.begin( )->first;
            const double finalTime = numericalIntegrationResult.rbegin( )->first;

            // Extract states
            const Eigen::Vector6d initialState = numericalIntegrationResult.begin( )->second;
            const Eigen::Vector6d finalState = numericalIntegrationResult.rbegin( )->second;

            // Create navigation system
            Eigen::Vector9d initialEstimatedStateVector = Eigen::Vector9d::Zero( );
            initialEstimatedStateVector.segment( NavigationSystem::cartesian_position_index, 6 ) = initialState;
            boost::shared_ptr< NavigationSystem > navigationSystem =
                    systemGenerator->createNavigationSystem( initialEstimatedStateVector, onboardTimeStep );
            const double standardGravitationalParameter = navigationSystem->getStandardGravitationalParameter( );
            const double planetaryRadius = navigationSystem->getRadius( );

            // Extract actual periapsis state
            unsigned int i = 0;
            Eigen::VectorXd historyOfActualAltitudes;
            Eigen::VectorXd historyOfActualVelocities;
            historyOfActualAltitudes.resize( numericalIntegrationResult.size( ) );
            historyOfActualVelocities.resizeLike( historyOfActualAltitudes );
            for ( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = numericalIntegrationResult.begin( );
                  stateIterator != numericalIntegrationResult.end( ); stateIterator++, i++ )
            {
                historyOfActualAltitudes[ i ] = stateIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
                historyOfActualVelocities[ i ] = stateIterator->second.segment( 3, 3 ).norm( );
            }
            Eigen::VectorXd::Index actualPeriapsisIndex;
            double actualPeriapsisAltitude = historyOfActualAltitudes.minCoeff( &actualPeriapsisIndex );
            double actualPeriapsisVelocity = historyOfActualVelocities[ actualPeriapsisIndex ];

            // Loop over various test cases
            for ( unsigned int navigationMode = testModes.first; navigationMode < testModes.second; navigationMode++ )
            {
                // Inform user
                std::cout << std::endl << "=========================================================================" << std::endl
                          << "Sun gravity: " << sunGravityIndex << ". Initial conditions: " << initialConditions
                          << ". Navigation mode: " << navigationMode << "." << std::endl
                          << "=========================================================================" << std::endl;

                // Determine initial state
                Eigen::Vector6d estimatedInitialState = initialState;
                switch ( navigationMode )
                {
                case 1:
                {
                    // Generate noise distributions
                    boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                            statistics::createBoostContinuousRandomVariableGenerator(
                                statistics::normal_boost_distribution, { 0.0, 50.0 }, 1 );
                    boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                            statistics::createBoostContinuousRandomVariableGenerator(
                                statistics::normal_boost_distribution, { 0.0, 0.1 }, 2 );

                    // Add random noise
                    for ( unsigned int i = 0; i < 6; i++ )
                    {
                        if ( i < 3 )
                        {
                            estimatedInitialState[ i ] += positionNoiseGenerator->getRandomVariableValue( );
                        }
                        else
                        {
                            estimatedInitialState[ i ] += velocityNoiseGenerator->getRandomVariableValue( );
                        }
                    }
                    break;
                }
                }

                // Set navigation conditions
                navigationSystem->setCurrentTime( initialTime );
                navigationSystem->setCurrentEstimatedCartesianState( estimatedInitialState );

                // Propagate state with onboard software
                boost::shared_ptr< propagators::PropagationTimeTerminationSettings > terminationSettings =
                        boost::make_shared< propagators::PropagationTimeTerminationSettings >( finalTime );
                const std::map< double, Eigen::VectorXd > estimatedResult =
                        navigationSystem->propagateTranslationalStateWithCustomTerminationSettings( terminationSettings ).second.first;

                // Extract estimated periapsis state
                unsigned int i = 0;
                Eigen::VectorXd historyOfEstimatedAltitudes;
                Eigen::VectorXd historyOfEstimatedVelocities;
                historyOfEstimatedAltitudes.resize( estimatedResult.size( ) );
                historyOfEstimatedVelocities.resizeLike( historyOfEstimatedAltitudes );
                for ( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = estimatedResult.begin( );
                      stateIterator != estimatedResult.end( ); stateIterator++, i++ )
                {
                    historyOfEstimatedAltitudes[ i ] = stateIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
                    historyOfEstimatedVelocities[ i ] = stateIterator->second.segment( 3, 3 ).norm( );
                }
                Eigen::VectorXd::Index estimatedPeriapsisIndex;
                double estimatedPeriapsisAltitude = historyOfEstimatedAltitudes.minCoeff( &estimatedPeriapsisIndex );
                double estimatedPeriapsisVelocity = historyOfEstimatedVelocities[ estimatedPeriapsisIndex ];

                // Store results
                propagationResults[ sunGravityIndex ][ initialConditions ][ navigationMode ] =
                        ( Eigen::VectorXd( 4 ) << actualPeriapsisAltitude, estimatedPeriapsisAltitude,
                          actualPeriapsisVelocity, estimatedPeriapsisVelocity ).finished( );

                // Compare results
                std::cout << "Actual: " << actualPeriapsisAltitude / 1.0e3 << " km. "
                          << "Estimated: " << estimatedPeriapsisAltitude / 1.0e3 << " km." << std::endl;

                // Show final state
                if ( ( sunGravityIndex == 0 ) && ( navigationMode == 0 ) )
                {
                    Eigen::VectorXd finalStateInKeplerianElements =
                            orbital_element_conversions::convertCartesianToKeplerianElements< double >( estimatedResult.rbegin( )->second,
                                                                                                        standardGravitationalParameter );
                    std::cout << std::setprecision( 16 ) << "Final state: "
                              << finalStateInKeplerianElements.transpose( ) << std::setprecision( 6 ) << std::endl;
                }
            }
        }
    }

    // Save results to file
    std::map< int, std::string > outFileNamesMap;
    outFileNamesMap[ 0 ] = "TestingResults/nsSeActualAltitude.dat";
    outFileNamesMap[ 1 ] = "TestingResults/nsSeEstimatedAltitude.dat";
    outFileNamesMap[ 2 ] = "TestingResults/nsSeActualVelocity.dat";
    outFileNamesMap[ 3 ] = "TestingResults/nsSeEstimatedVelocity.dat";
    std::vector< double > vectorOfSuns = { 0.0, 1.0 };
    Eigen::VectorXd vectorOfConditions = Eigen::ArrayXd::LinSpaced( ( testConditions.second - testConditions.first ),
                                                                    testConditions.first, testConditions.second - 1 );
    Eigen::VectorXd vectorOfModes = Eigen::ArrayXd::LinSpaced( ( testModes.second - testModes.first ),
                                                               testModes.first, testModes.second - 1 );
    std::vector< std::vector< double > > independentVariables;
    independentVariables.push_back( vectorOfSuns );
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfConditions ) );
    independentVariables.push_back( utilities::convertEigenVectorToStlVector( vectorOfModes ) );
    input_output::MultiArrayFileWriter< 3, 4 >::writeMultiArrayAndIndependentVariablesToFiles( outFileNamesMap,
                                                                                               independentVariables,
                                                                                               propagationResults );
}

//BOOST_AUTO_TEST_CASE( testAccelerometerPostProcessing )
//{

//}

BOOST_AUTO_TEST_CASE( testPeriapseTimeEstimator )
{
    // Declare navigation system generator
    boost::shared_ptr< NavigationSystemGenerator > systemGenerator =
            boost::make_shared< NavigationSystemGenerator >( aerodynamics::exponential_atmosphere_model );

    // Define constants
    const double onboardTimeStep = 0.1;

    // Loop over various initial conditions
    for ( unsigned int initialConditions = testConditions.first; initialConditions < testConditions.second; initialConditions++ )
    {
        // Propagate spacecraft state and retireve solution
        const std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > fullNumericalResults =
                propagateStateForAerobrakingScenario( initialConditions );
        const std::map< double, Eigen::VectorXd > numericalIntegrationResult = fullNumericalResults.first;
        const std::map< double, Eigen::VectorXd > dependentVariablesResult = fullNumericalResults.second;

        // Extract times
        const double initialTime = numericalIntegrationResult.begin( )->first;
        const double finalTime = numericalIntegrationResult.rbegin( )->first;

        // Extract states
        const Eigen::Vector6d initialState = numericalIntegrationResult.begin( )->second;
        const Eigen::Vector6d finalState = numericalIntegrationResult.rbegin( )->second;

        // Create navigation system
        Eigen::Vector9d initialEstimatedStateVector = Eigen::Vector9d::Zero( );
        initialEstimatedStateVector.segment( NavigationSystem::cartesian_position_index, 6 ) = initialState;
        boost::shared_ptr< NavigationSystem > navigationSystem =
                systemGenerator->createNavigationSystem( initialEstimatedStateVector, onboardTimeStep );
        const double standardGravitationalParameter = navigationSystem->getStandardGravitationalParameter( );
        const double planetaryRadius = navigationSystem->getRadius( );

        // Extract actual periapsis state
        unsigned int i = 0;
        Eigen::VectorXd historyOfActualAltitudes;
        Eigen::VectorXd historyOfTimes;
        historyOfActualAltitudes.resize( numericalIntegrationResult.size( ) );
        historyOfTimes.resizeLike( historyOfActualAltitudes );
        for ( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = numericalIntegrationResult.begin( );
              stateIterator != numericalIntegrationResult.end( ); stateIterator++, i++ )
        {
            historyOfActualAltitudes[ i ] = stateIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
            historyOfTimes[ i ] = stateIterator->first;
        }
        Eigen::VectorXd::Index actualPeriapsisIndex;
        double actualPeriapsisAltitude = historyOfActualAltitudes.minCoeff( &actualPeriapsisIndex );
        double actualPeriapsisTime = historyOfTimes[ actualPeriapsisIndex ];

        // Interpolate to get result at constant time step
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > numericalIntegrationInterpolator =
                boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >( numericalIntegrationResult, 8 );
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariableInterpolator =
                boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >( dependentVariablesResult, 8 );

        // Loop over various test cases
        for ( unsigned int navigationMode = testModes.first; navigationMode < testModes.second; navigationMode++ )
        {
            // Inform user
            std::cout << std::endl << "========================================================" << std::endl
                      << "Initial conditions: " << initialConditions << ". Test case: " << navigationMode << "." << std::endl
                      << "========================================================" << std::endl;

            // Accelerometer noise
            double accelerometerAccuracy = 2.0e-4 * std::sqrt( 1 / 0.02 ) * 0.1;
            boost::shared_ptr< statistics::RandomVariableGenerator< double > > accelerationNoiseGenerator =
                    statistics::createBoostContinuousRandomVariableGenerator(
                        statistics::normal_boost_distribution, { 0.0, accelerometerAccuracy }, 0 );

            // Get estimated state based on test case
            double timeAtDAIA;
            std::vector< double > aerodynamicAcceleration;
            std::map< double, Eigen::VectorXd > estimatedResults;
            std::map< double, Eigen::VectorXd > dependentVariables;
            switch ( navigationMode )
            {
            case 0:
            {
                // Generate noise distributions
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 50.0 }, 1 );
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 0.1 }, 2 );

                // Loop over each numerical result
                double currentTime = initialTime;
                while ( currentTime < finalTime )
                {
                    // Interpolate at current time
                    estimatedResults[ currentTime ] = numericalIntegrationInterpolator->interpolate( currentTime );
                    dependentVariables[ currentTime ] = dependentVariableInterpolator->interpolate( currentTime );

                    // Add random noise
//                    for ( unsigned int i = 0; i < 6; i++ )
//                    {
//                        if ( i < 3 )
//                        {
//                            estimatedResults[ currentTime ][ i ] += positionNoiseGenerator->getRandomVariableValue( );
//                        }
//                        else
//                        {
//                            estimatedResults[ currentTime ][ i ] += velocityNoiseGenerator->getRandomVariableValue( );
//                        }
//                    }

                    // Set state in navigation system
                    navigationSystem->setCurrentTime( currentTime );
                    navigationSystem->setCurrentEstimatedCartesianState( estimatedResults[ currentTime ] );

                    // Compute aerodynamic acceleration norm
                    Eigen::Vector3d currentAerodynamicAcceleration = dependentVariables[ currentTime ].segment( 0, 3 );
                    for ( unsigned int i = 0; i < 3; i++ )
                    {
                        currentAerodynamicAcceleration[ i ] += accelerationNoiseGenerator->getRandomVariableValue( );
                    }
                    aerodynamicAcceleration.push_back( currentAerodynamicAcceleration.norm( ) );

                    // Check whether it is time to stop propagation and run post-atmosphere processes
                    double currentEstimatedTrueAnomaly = navigationSystem->getCurrentEstimatedTranslationalState( ).second[ 5 ];
                    if ( navigationSystem->getIsSpacecraftAboveDynamicAtmosphericInterfaceAltitude( ) &&
                         ( currentEstimatedTrueAnomaly < ( 0.95 * mathematical_constants::PI ) ) )
                    {
                        timeAtDAIA = currentTime;
                        break;
                    }

                    // Step time forward
                    currentTime += onboardTimeStep;
                }
                break;
            }
            case 1:
            {
                // Generate noise distributions
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 0.5 }, 1 );
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 0.001 }, 2 );
                boost::function< double( const double, const double ) > offsetFunction =
                        boost::bind( &offsetFunctionWrapper, initialTime, _1, _2 );

                // Loop over each numerical result
                double currentTime = initialTime;
                while ( currentTime < finalTime )
                {
                    // Interpolate at current time
                    estimatedResults[ currentTime ] = numericalIntegrationInterpolator->interpolate( currentTime );
                    dependentVariables[ currentTime ] = dependentVariableInterpolator->interpolate( currentTime );

                    // Compute aerodynamic acceleration norm
                    Eigen::Vector3d currentAerodynamicAcceleration = dependentVariables[ currentTime ].segment( 0, 3 );
                    for ( unsigned int i = 0; i < 3; i++ )
                    {
                        currentAerodynamicAcceleration[ i ] += accelerationNoiseGenerator->getRandomVariableValue( );
                    }
                    aerodynamicAcceleration.push_back( currentAerodynamicAcceleration.norm( ) );

                    // Add random noise
                    double currentOffset;
                    for ( unsigned int i = 0; i < 6; i++ )
                    {
                        if ( i < 3 )
                        {
                            currentOffset = offsetFunction( currentTime, dependentVariables[ currentTime ][ i ] );
                            estimatedResults[ currentTime ][ i ] += currentOffset + positionNoiseGenerator->getRandomVariableValue( );
                        }
                        else
                        {
                            estimatedResults[ currentTime ][ i ] += velocityNoiseGenerator->getRandomVariableValue( );
                        }
                    }

                    // Set state in navigation system
                    navigationSystem->setCurrentTime( currentTime );
                    navigationSystem->setCurrentEstimatedCartesianState( estimatedResults[ currentTime ] );

                    // Check whether it is time to stop propagation and run post-atmosphere processes
                    double currentEstimatedTrueAnomaly = navigationSystem->getCurrentEstimatedTranslationalState( ).second[ 5 ];
                    if ( navigationSystem->getIsSpacecraftAboveDynamicAtmosphericInterfaceAltitude( ) &&
                         ( currentEstimatedTrueAnomaly < ( 0.95 * mathematical_constants::PI ) ) )
                    {
                        timeAtDAIA = currentTime;
                        break;
                    }

                    // Step time forward
                    currentTime += onboardTimeStep;
                }
                break;
            }
            default:
                throw std::runtime_error( "This mode is not supported yet." );
            }

            // Extract data below atmospheric interface
            unsigned int i = 0;
            Eigen::VectorXd historyOfAltitudes;
            std::map< double, Eigen::Vector6d > mapOfEstimatedKeplerianStatesBelowAtmosphericInterface;
            std::vector< double > vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface;
            for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = estimatedResults.begin( );
                  mapIterator != estimatedResults.end( ); mapIterator++, i++ )
            {
                if ( mapIterator->second.segment( 0, 3 ).norm( ) <= navigationSystem->getAtmosphericInterfaceRadius( ) )
                {
                    mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ mapIterator->first ] =
                            orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                                mapIterator->second, standardGravitationalParameter );
                    vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back( aerodynamicAcceleration.at( i ) );
                }
            }

            // Run periapse time estimator based on current data
            Eigen::Vector6d estimatedChangeInKeplerianElements =
                    navigationSystem->testPeriapseTimeEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                                                 vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

            // Compute actual change in elements
            Eigen::Vector6d actualChangeInKeplerianElements =
                    orbital_element_conversions::convertCartesianToKeplerianElements( finalState, standardGravitationalParameter ) -
                    orbital_element_conversions::convertCartesianToKeplerianElements( initialState, standardGravitationalParameter );

            // Check how close the estimates are
            std::cout << "Actual: ";
            const std::vector< unsigned int > indicesToCheck = { 0, 1 };
            for ( unsigned int i: indicesToCheck )
            {
                std::cout << actualChangeInKeplerianElements[ i ] << " ";
                BOOST_CHECK_CLOSE_FRACTION( estimatedChangeInKeplerianElements[ i ], actualChangeInKeplerianElements[ i ], 1.0e-1 );
            }
            std::cout << std::endl << "Actual periapsis time: " << actualPeriapsisTime - initialTime << " s" << std::endl;

            // Check true anomaly estimate by looking at shift in orbit
            Eigen::Vector6d propagatedKeplerianStateAtDAIA = orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                        numericalIntegrationInterpolator->interpolate( timeAtDAIA ), standardGravitationalParameter );
            Eigen::Vector6d estimatedKeplerianStateAtDAIA = orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                        estimatedResults.rbegin( )->second, standardGravitationalParameter );
            BOOST_CHECK_CLOSE_FRACTION( estimatedKeplerianStateAtDAIA[ 5 ] - propagatedKeplerianStateAtDAIA[ 5 ],
                    estimatedChangeInKeplerianElements[ 5 ], 1.0e-1 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testAtmosphereEstimator )
{
    // Loop over various atmosphere models
    for ( unsigned int atmosphereModel = testAtmospheres.first; atmosphereModel < testAtmospheres.second; atmosphereModel++ )
    {
        if ( atmosphereModel != 1 )
        {
            // Declare navigation system generator
            boost::shared_ptr< NavigationSystemGenerator > systemGenerator =
                    boost::make_shared< NavigationSystemGenerator >(
                        static_cast< aerodynamics::AvailableConstantTemperatureAtmosphereModels >( atmosphereModel ) );

            // Define constants
            const double onboardTimeStep = 0.1;

            // Loop over various initial conditions
            for ( unsigned int initialConditions = testConditions.first; initialConditions < testConditions.second; initialConditions++ )
            {
                // Propagate spacecraft state and retireve solution
                const std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > fullNumericalResults =
                        propagateStateForAerobrakingScenario( initialConditions );
                const std::map< double, Eigen::VectorXd > numericalIntegrationResult = fullNumericalResults.first;
                const std::map< double, Eigen::VectorXd > dependentVariablesResult = fullNumericalResults.second;

                // Extract times
                const double initialTime = numericalIntegrationResult.begin( )->first;
                const double finalTime = numericalIntegrationResult.rbegin( )->first;

                // Extract states
                const Eigen::Vector6d initialState = numericalIntegrationResult.begin( )->second;
                const Eigen::Vector6d finalState = numericalIntegrationResult.rbegin( )->second;

                // Create navigation system
                Eigen::Vector9d initialEstimatedStateVector = Eigen::Vector9d::Zero( );
                initialEstimatedStateVector.segment( NavigationSystem::cartesian_position_index, 6 ) = initialState;
                boost::shared_ptr< NavigationSystem > navigationSystem =
                        systemGenerator->createNavigationSystem( initialEstimatedStateVector, onboardTimeStep );
                const double standardGravitationalParameter = navigationSystem->getStandardGravitationalParameter( );

                // Interpolate to get result at constant time step
                boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > numericalIntegrationInterpolator =
                        boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >(
                            numericalIntegrationResult, 8 );
                boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariableInterpolator =
                        boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >(
                            dependentVariablesResult, 8 );

                // Loop over various test cases
                for ( unsigned int navigationMode = testModes.first; navigationMode < testModes.second; navigationMode++ )
                {
                    // Inform user
                    std::cout << std::endl << "=========================================================================" << std::endl
                              << "Atmoshere Model: " << atmosphereModel << ". Initial conditions: " << initialConditions
                              << ". Test case: " << navigationMode << "." << std::endl
                              << "=========================================================================" << std::endl;

                    // Accelerometer noise
                    double accelerometerAccuracy = 2.0e-4 * std::sqrt( 1 / 0.02 ) * 0.1;
                    boost::shared_ptr< statistics::RandomVariableGenerator< double > > accelerationNoiseGenerator =
                            statistics::createBoostContinuousRandomVariableGenerator(
                                statistics::normal_boost_distribution, { 0.0, accelerometerAccuracy }, 0 );

                    // Get estimated state based on test case
                    double timeAtDAIA;
                    std::vector< double > aerodynamicAcceleration;
                    std::map< double, Eigen::VectorXd > numericalResults;
                    std::map< double, Eigen::VectorXd > estimatedResults;
                    std::map< double, Eigen::VectorXd > dependentVariables;
                    switch ( navigationMode )
                    {
                    case 0:
                    {
                        // Generate noise distributions
                        boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                                statistics::createBoostContinuousRandomVariableGenerator(
                                    statistics::normal_boost_distribution, { 0.0, 50.0 }, 1 );
                        boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                                statistics::createBoostContinuousRandomVariableGenerator(
                                    statistics::normal_boost_distribution, { 0.0, 0.1 }, 2 );

                        // Loop over each numerical result
                        double currentTime = initialTime;
                        while ( currentTime < finalTime )
                        {
                            // Interpolate at current time
                            numericalResults[ currentTime ] = numericalIntegrationInterpolator->interpolate( currentTime );
                            estimatedResults[ currentTime ] = numericalResults[ currentTime ];
                            dependentVariables[ currentTime ] = dependentVariableInterpolator->interpolate( currentTime );

                            // Add random noise
                            for ( unsigned int i = 0; i < 6; i++ )
                            {
                                if ( i < 3 )
                                {
                                    estimatedResults[ currentTime ][ i ] += positionNoiseGenerator->getRandomVariableValue( );
                                }
                                else
                                {
                                    estimatedResults[ currentTime ][ i ] += velocityNoiseGenerator->getRandomVariableValue( );
                                }
                            }

                            // Set state in navigation system
                            navigationSystem->setCurrentTime( currentTime );
                            navigationSystem->setCurrentEstimatedCartesianState( estimatedResults[ currentTime ] );

                            // Compute aerodynamic acceleration norm
                            Eigen::Vector3d currentAerodynamicAcceleration = dependentVariables[ currentTime ].segment( 0, 3 );
                            for ( unsigned int i = 0; i < 3; i++ )
                            {
                                currentAerodynamicAcceleration[ i ] += accelerationNoiseGenerator->getRandomVariableValue( );
                            }
                            aerodynamicAcceleration.push_back( currentAerodynamicAcceleration.norm( ) );

                            // Check whether it is time to stop propagation and run post-atmosphere processes
                            double currentEstimatedTrueAnomaly = navigationSystem->getCurrentEstimatedTranslationalState( ).second[ 5 ];
                            if ( navigationSystem->getIsSpacecraftAboveDynamicAtmosphericInterfaceAltitude( ) &&
                                 ( currentEstimatedTrueAnomaly < ( 0.95 * mathematical_constants::PI ) ) )
                            {
                                timeAtDAIA = currentTime;
                                break;
                            }

                            // Step time forward
                            currentTime += onboardTimeStep;
                        }
                        break;
                    }
                    case 1:
                    {
                        // Generate noise distributions
                        boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                                statistics::createBoostContinuousRandomVariableGenerator(
                                    statistics::normal_boost_distribution, { 0.0, 0.5 }, 1 );
                        boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                                statistics::createBoostContinuousRandomVariableGenerator(
                                    statistics::normal_boost_distribution, { 0.0, 0.001 }, 2 );
                        boost::function< double( const double, const double ) > offsetFunction =
                                boost::bind( &offsetFunctionWrapper, initialTime, _1, _2 );

                        // Loop over each numerical result
                        double currentTime = initialTime;
                        while ( currentTime < finalTime )
                        {
                            // Interpolate at current time
                            numericalResults[ currentTime ] = numericalIntegrationInterpolator->interpolate( currentTime );
                            estimatedResults[ currentTime ] = numericalResults[ currentTime ];
                            dependentVariables[ currentTime ] = dependentVariableInterpolator->interpolate( currentTime );

                            // Compute aerodynamic acceleration norm
                            Eigen::Vector3d currentAerodynamicAcceleration = dependentVariables[ currentTime ].segment( 0, 3 );
                            for ( unsigned int i = 0; i < 3; i++ )
                            {
                                currentAerodynamicAcceleration[ i ] += accelerationNoiseGenerator->getRandomVariableValue( );
                            }
                            aerodynamicAcceleration.push_back( currentAerodynamicAcceleration.norm( ) );

                            // Add random noise
                            double currentOffset;
                            for ( unsigned int i = 0; i < 6; i++ )
                            {
                                if ( i < 3 )
                                {
                                    currentOffset = offsetFunction( currentTime, dependentVariables[ currentTime ][ i ] );
                                    estimatedResults[ currentTime ][ i ] += currentOffset + positionNoiseGenerator->getRandomVariableValue( );
                                }
                                else
                                {
                                    estimatedResults[ currentTime ][ i ] += velocityNoiseGenerator->getRandomVariableValue( );
                                }
                            }

                            // Set state in navigation system
                            navigationSystem->setCurrentTime( currentTime );
                            navigationSystem->setCurrentEstimatedCartesianState( estimatedResults[ currentTime ] );

                            // Check whether it is time to stop propagation and run post-atmosphere processes
                            double currentEstimatedTrueAnomaly = navigationSystem->getCurrentEstimatedTranslationalState( ).second[ 5 ];
                            if ( navigationSystem->getIsSpacecraftAboveDynamicAtmosphericInterfaceAltitude( ) &&
                                 ( currentEstimatedTrueAnomaly < ( 0.95 * mathematical_constants::PI ) ) )
                            {
                                timeAtDAIA = currentTime;
                                break;
                            }

                            // Step time forward
                            currentTime += onboardTimeStep;
                        }
                        break;
                    }
                    default:
                        throw std::runtime_error( "This mode is not supported yet." );
                    }

                    // Extract data below atmospheric interface
                    unsigned int i = 0;
                    std::map< double, Eigen::Vector6d > mapOfActualDensityBelowAtmosphericInterface;
                    std::map< double, Eigen::Vector6d > mapOfActualKeplerianStatesBelowAtmosphericInterface;
                    std::map< double, Eigen::Vector6d > mapOfEstimatedKeplerianStatesBelowAtmosphericInterface;
                    std::vector< double > vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface;
                    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = estimatedResults.begin( );
                          mapIterator != estimatedResults.end( ); mapIterator++, i++ )
                    {
                        if ( mapIterator->second.segment( 0, 3 ).norm( ) <= navigationSystem->getAtmosphericInterfaceRadius( ) )
                        {
                            mapOfActualKeplerianStatesBelowAtmosphericInterface[ mapIterator->first ] =
                                    orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                                        numericalResults[ mapIterator->first ], standardGravitationalParameter );
                            mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ mapIterator->first ] =
                                    orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                                        mapIterator->second, standardGravitationalParameter );
                            vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back( aerodynamicAcceleration.at( i ) );
                        }
                    }
                    input_output::writeDataMapToTextFile( mapOfActualKeplerianStatesBelowAtmosphericInterface,
                                                          "nsKepler_act.dat", "TestingResults/" );

                    // Run periapse time estimator based on current data
                    Eigen::VectorXd estimatedAtmosphereParameters =
                            navigationSystem->testAtmosphereEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                                                       vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

                    // Create onboard atmosphere model with estimated parameters
                    boost::shared_ptr< aerodynamics::CustomConstantTemperatureAtmosphere > onboardAtmosphereModel =
                            boost::make_shared< aerodynamics::CustomConstantTemperatureAtmosphere >(
                                static_cast< aerodynamics::AvailableConstantTemperatureAtmosphereModels >( atmosphereModel ),
                                215.0, 197.0, 1.3, utilities::convertEigenVectorToStlVector( estimatedAtmosphereParameters ) );

                    // Loop over conditions and extract density
                    i = 0;
                    double currentDensity;
                    double estimatedDensity;
                    double currentRadialDistance;
                    std::vector< double > errorValues;
                    std::vector< double > altitudeValues;
                    std::map< double, Eigen::Vector3d > densityValues;
                    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = numericalResults.begin( );
                          mapIterator != numericalResults.end( ); mapIterator++, i++ )
                    {
                        currentRadialDistance = mapIterator->second.segment( 0, 3 ).norm( );
                        if ( currentRadialDistance <= navigationSystem->getReducedAtmosphericInterfaceRadius( ) )
                        {
                            currentDensity = dependentVariables[ mapIterator->first ][ 3 ];
                            estimatedDensity = onboardAtmosphereModel->getDensity( currentRadialDistance - navigationSystem->getRadius( ) );
                            densityValues[ mapIterator->first ] =
                                    ( Eigen::VectorXd( 3 ) << currentRadialDistance - navigationSystem->getRadius( ),
                                      currentDensity, estimatedDensity ).finished( );
                            errorValues.push_back( ( estimatedDensity - currentDensity ) / currentDensity );
                            altitudeValues.push_back( currentRadialDistance - navigationSystem->getRadius( ) );
                        }
                    }
                    if ( ( ( atmosphereModel == 0 ) && ( initialConditions == 2 ) && ( navigationMode == 0 ) ) ||
                         ( ( atmosphereModel == 2 ) && ( initialConditions == 2 ) && ( navigationMode == 0 ) ) )
                    {
                        input_output::writeDataMapToTextFile( densityValues, "nsAeDensity_" + std::to_string( atmosphereModel ) + ".dat",
                                                              "TestingResults/" );
                    }

                    Eigen::VectorXd::Index indexMaximum;
                    Eigen::VectorXd::Index indexMinumum;
                    double maximumError = utilities::convertStlVectorToEigenVector( errorValues ).maxCoeff( &indexMaximum );
                    double minimumError = utilities::convertStlVectorToEigenVector( errorValues ).minCoeff( &indexMinumum );
                    Eigen::VectorXd altitudeValuesX = utilities::convertStlVectorToEigenVector( altitudeValues );

                    // Check that resulting atmospheres are similar around the periapsis area
                    BOOST_CHECK_LE( maximumError, 0.15 );
                    BOOST_CHECK_GE( minimumError, -0.15 );
                    std::cout << "Max. error: " << maximumError << " at " << altitudeValuesX[ indexMaximum ] / 1.0e3 << " km" << std::endl
                              << "Min. error: " << minimumError << " at " << altitudeValuesX[ indexMinumum ] / 1.0e3 << " km" << std::endl;
                }
            }
        }
    }
}

//BOOST_AUTO_TEST_CASE( testDSNTracking )
//{

//}

BOOST_AUTO_TEST_SUITE_END( )

}

}

