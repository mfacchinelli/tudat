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
    for ( unsigned int initialConditions = 0; initialConditions < 1; initialConditions++ )
    {
        // Propagate spacecraft state and retireve solution
        const std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > fullNumericalResults =
                propagateStateForAerobrakingScenario( initialConditions );
        const std::map< double, Eigen::VectorXd > numericalIntegrationResult = fullNumericalResults.first;
        const std::map< double, Eigen::VectorXd > dependentVariablesResult = fullNumericalResults.second;

        // Extract times
        const double initialTime = numericalIntegrationResult.begin( )->first;
        const double finalTime = numericalIntegrationResult.rbegin( )->first;

        // Create navigation system
        Eigen::Vector9d initialEstimatedStateVector = Eigen::Vector9d::Zero( );
        initialEstimatedStateVector.segment( NavigationSystem::cartesian_position_index, 6 ) = numericalIntegrationResult.begin( )->second;
        boost::shared_ptr< NavigationSystem > navigationSystem =
                systemGenerator->createNavigationSystem( initialEstimatedStateVector, onboardTimeStep );

        // Interpolate to get result at constant time step
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > numericalIntegrationInterpolator =
                boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >( numericalIntegrationResult, 8 );
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariableInterpolator =
                boost::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd, double > >( dependentVariablesResult, 8 );

        // Loop over various test cases
        for ( unsigned int testCase = 0; testCase < 1; testCase++ )
        {
            // Get estimated state based on test case
            std::vector< double > aerodynamicAcceleration;
            std::map< double, Eigen::VectorXd > estimatedResults;
            std::map< double, Eigen::VectorXd > dependentVariables;
            switch ( testCase )
            {
            case 0:
            {
                // Generate noise distributions
                double accelerometerAccuracy = 2.0e-4 * std::sqrt( 1 / 0.02 ) * 0.1;
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > positionNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 500.0 }, 1 );
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > velocityNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, 10.0 }, 1 );
                boost::shared_ptr< statistics::RandomVariableGenerator< double > > accelerationNoiseGenerator =
                        statistics::createBoostContinuousRandomVariableGenerator(
                            statistics::normal_boost_distribution, { 0.0, accelerometerAccuracy }, 1 );

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

                    // Compute aerodynamic acceleration norm
                    Eigen::Vector3d currentAerodynamicAcceleration = dependentVariables[ currentTime ].segment( 0, 3 );
                    for ( unsigned int i = 0; i < 3; i++ )
                    {
                        currentAerodynamicAcceleration[ i ] += accelerationNoiseGenerator->getRandomVariableValue( );
                    }
                    aerodynamicAcceleration.push_back( currentAerodynamicAcceleration.norm( ) );

                    // Step time forward
                    currentTime += onboardTimeStep;
                }
                break;
            }
            case 1:
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
                                mapIterator->second, navigationSystem->getStandardGravitationalParameter( ) );
                    vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back( aerodynamicAcceleration.at( i ) );
                }
            }

            // Run periapse time estimator based on current data
            navigationSystem->testPeriapseTimeEstimator( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface,
                                                         vectorOfMeasuredAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
        }
    }
}


//BOOST_AUTO_TEST_CASE( testAtmosphereEstimator )
//{

//}

//BOOST_AUTO_TEST_CASE( testDSNTracking )
//{

//}

BOOST_AUTO_TEST_SUITE_END( )

}

}

