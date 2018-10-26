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
const std::pair< unsigned int, unsigned int > testConditions = { 0, 2 };
const std::pair< unsigned int, unsigned int > testModes = { 0, 2 };

BOOST_AUTO_TEST_SUITE( test_navigation_system )

//BOOST_AUTO_TEST_CASE( testNavigationPhase )
//{

//}

//BOOST_AUTO_TEST_CASE( testStateEstimator )
//{

//}

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
        const double standardGravitationalParameters = navigationSystem->getStandardGravitationalParameter( );

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
                            statistics::normal_boost_distribution, { 0.0, 1.0 }, 2 );

                // Loop over each numerical result
                double currentTime = initialTime;
                while ( currentTime < finalTime )
                {
                    // Interpolate at current time
                    estimatedResults[ currentTime ] = numericalIntegrationInterpolator->interpolate( currentTime );
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
                            statistics::normal_boost_distribution, { 0.0, 5.0 }, 1 );
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 0.1 }, 2 );

                // Linear offset function
                std::function< double( const double, const double ) > offsetFunction =
                        [ = ]( const double time, const double acceleration ){
                    return ( 1.0e-2 + 2.5e1 * acceleration ) * ( time - initialTime ); };

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
            std::map< double, Eigen::Vector6d > mapOfEstimatedKeplerianStatesBelowAtmosphericInterface;
            std::vector< double > vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface;
            for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = estimatedResults.begin( );
                  mapIterator != estimatedResults.end( ); mapIterator++, i++ )
            {
                if ( mapIterator->second.segment( 0, 3 ).norm( ) <= navigationSystem->getAtmosphericInterfaceRadius( ) )
                {
                    mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ mapIterator->first ] =
                            orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                                mapIterator->second, standardGravitationalParameters );
                    vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back( aerodynamicAcceleration.at( i ) );
                }
            }

            // Run periapse time estimator based on current data
            Eigen::Vector6d estimatedChangeInKeplerianElements =
                    navigationSystem->testPeriapseTimeEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                                                 vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

            // Compute actual change in elements
            Eigen::Vector6d actualChangeInKeplerianElements =
                    orbital_element_conversions::convertCartesianToKeplerianElements( finalState, standardGravitationalParameters ) -
                    orbital_element_conversions::convertCartesianToKeplerianElements( initialState, standardGravitationalParameters );

            std::cout << "Previous: " << orbital_element_conversions::convertCartesianToKeplerianElements(
                             initialState, standardGravitationalParameters ).transpose( ) << std::endl
                      << "Current: " << orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                             numericalIntegrationInterpolator->interpolate( timeAtDAIA ),
                             standardGravitationalParameters ).transpose( ) << std::endl
                      << "Next: " << orbital_element_conversions::convertCartesianToKeplerianElements(
                             finalState, standardGravitationalParameters ).transpose( ) << std::endl;

            // Check how close the estimates are
            std::cout << "Actual: ";
            const std::vector< unsigned int > indicesToCheck = { 0, 1 };
            for ( unsigned int i: indicesToCheck )
            {
                std::cout << actualChangeInKeplerianElements[ i ] << " ";
                BOOST_CHECK_CLOSE_FRACTION( estimatedChangeInKeplerianElements[ i ], actualChangeInKeplerianElements[ i ], 1.0e-1 );
            }
            std::cout << std::endl;

            // Check true anomaly estimate by looking at shift in orbit
            Eigen::Vector6d propagatedKeplerianStateAtDAIA = orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                        numericalIntegrationInterpolator->interpolate( timeAtDAIA ), standardGravitationalParameters );
            Eigen::Vector6d estimatedKeplerianStateAtDAIA = orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                        estimatedResults.rbegin( )->second, standardGravitationalParameters );
            BOOST_CHECK_CLOSE_FRACTION( estimatedKeplerianStateAtDAIA[ 5 ] - propagatedKeplerianStateAtDAIA[ 5 ],
                    estimatedChangeInKeplerianElements[ 5 ], 1.0e-1 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testAtmosphereEstimator )
{
    // Expected results
    std::map< std::tuple< unsigned int, unsigned int, unsigned int >, std::pair< double, double > > expectedResults;
    expectedResults[ { 0, 0, 0 } ] = { -0.765026547209958, 1.31583846986147 };
    expectedResults[ { 0, 0, 1 } ] = { -0.740769780829461, 3.98142882969108 };
    expectedResults[ { 0, 1, 0 } ] = { -0.7548723576166, 2.54764862007878 };
    expectedResults[ { 0, 1, 1 } ] = { -0.799875652105309, 3.0534113745639 };
    expectedResults[ { 2, 0, 0 } ] = { -0.761742398031845, 1.31617390007452 };
    expectedResults[ { 2, 0, 1 } ] = { -0.731284595857299, 4.01316556176618 };
    expectedResults[ { 2, 1, 0 } ] = { -0.737521490108276, 2.72903166215168 };
    expectedResults[ { 2, 1, 1 } ] = { -0.784575862911139, 3.17126418088939 };

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
                const double standardGravitationalParameters = navigationSystem->getStandardGravitationalParameter( );

                // Interpolate to get result at constant time step
                boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > numericalIntegrationInterpolator =
                        boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >( numericalIntegrationResult, 8 );
                boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariableInterpolator =
                        boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >( dependentVariablesResult, 8 );

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
                                    statistics::normal_boost_distribution, { 0.0, 1.0 }, 2 );

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
                                    statistics::normal_boost_distribution, { 0.0, 5.0 }, 1 );
                        boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                                statistics::createBoostContinuousRandomVariableGenerator(
                                    statistics::normal_boost_distribution, { 0.0, 0.1 }, 2 );

                        // Linear offset function
                        std::function< double( const double, const double ) > offsetFunction =
                                [ = ]( const double time, const double acceleration ){
                            return ( 1.0e-2 + 2.5e1 * acceleration ) * ( time - initialTime ); };

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
                                        numericalResults[ mapIterator->first ], standardGravitationalParameters );
                            mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ mapIterator->first ] =
                                    orbital_element_conversions::convertCartesianToKeplerianElements< double >(
                                        mapIterator->second, standardGravitationalParameters );
                            vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back( aerodynamicAcceleration.at( i ) );
                        }
                    }
                    input_output::writeDataMapToTextFile( mapOfActualKeplerianStatesBelowAtmosphericInterface, "kepler_act.dat", "PTE&AEResults/" );

                    // Run periapse time estimator based on current data
                    Eigen::VectorXd estimatedAtmosphereParameters =
                            navigationSystem->testAtmosphereEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                                                       vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );

                    // Create onboard atmosphere model with estimated parameters
                    boost::shared_ptr< aerodynamics::CustomConstantTemperatureAtmosphere > onboardAtmosphereModel =
                            boost::make_shared< aerodynamics::CustomConstantTemperatureAtmosphere >(
                                static_cast< aerodynamics::AvailableConstantTemperatureAtmosphereModels >( atmosphereModel ),
                                215.0, 197.0, 1.3, utilities::convertEigenVectorToStlVector( estimatedAtmosphereParameters ) );

                    // Check that resulting atmospheres are similar around the periapsis area
                    i = 0;
                    double currentRadialDistance;
                    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = numericalResults.begin( );
                          mapIterator != numericalResults.end( ); mapIterator++, i++ )
                    {
                        currentRadialDistance = mapIterator->second.segment( 0, 3 ).norm( );
                        if ( currentRadialDistance <= navigationSystem->getReducedAtmosphericInterfaceRadius( ) )
                        {
                            BOOST_CHECK_CLOSE_FRACTION( dependentVariables[ mapIterator->first ][ 3 ],
                                    onboardAtmosphereModel->getDensity( currentRadialDistance - navigationSystem->getRadius( ) ), 0.3 );
                        }
                    }
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

