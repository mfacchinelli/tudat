/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Mathematics/Statistics/basicStatistics.h"
#include "Tudat/Mathematics/Filters/extendedKalmanFilter.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

#include "Tudat/Mathematics/Filters/UnitTests/controlClass.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_extended_kalman_filter )

// Functions for extended Kalman filter.
Eigen::Vector2d stateFunction1( const double time, const Eigen::Vector2d& state, const Eigen::Vector2d& control )
{
    Eigen::Vector2d stateDerivative;
    stateDerivative[ 0 ] = state[ 1 ] * std::pow( std::cos( state[ 0 ] ), 3 );
    stateDerivative[ 1 ] = std::sin( state[ 0 ] );
    return stateDerivative;
}
Eigen::Vector1d measurementFunction1( const double time, const Eigen::Vector2d& state )
{
    Eigen::Vector1d measurement;
    measurement[ 0 ] = std::pow( state[ 0 ], 3 );
    return measurement;
}
Eigen::Matrix2d stateJacobianFunction1( const double time, const Eigen::Vector2d& state, const Eigen::Vector2d& control )
{
    Eigen::Matrix2d stateJacobian = Eigen::Matrix2d::Zero( );
    stateJacobian( 0, 0 ) = - 3 * state[ 1 ] * std::pow( std::cos( state[ 0 ] ), 2 ) * std::sin( state[ 0 ] );
    stateJacobian( 0, 1 ) = std::pow( std::cos( state[ 0 ] ), 3 );
    stateJacobian( 1, 0 ) = std::cos( state[ 0 ] );
    return stateJacobian;
}
Eigen::RowVector2d measurementJacobianFunction1( const double time, const Eigen::Vector2d& state )
{
    Eigen::RowVector2d measurementJacobian;
    measurementJacobian[ 0 ] = 3.0 * std::pow( state[ 0 ], 2 );
    return measurementJacobian;
}

// Test implementation of extended Kalman filter class.
BOOST_AUTO_TEST_CASE( testExtendedKalmanFilterFirstCase )
{
    using namespace tudat::filters;

    // Set initial conditions
    const double initialTime = 0;
    const double timeStep = 0.01;
    const unsigned int numberOfTimeSteps = 1000;

    Eigen::Vector2d initialStateVector;
    initialStateVector[ 0 ] = 3.0;
    initialStateVector[ 1 ] = -0.3;

    Eigen::Vector2d initialEstimatedStateVector;
    initialEstimatedStateVector[ 0 ] = 10.0;
    initialEstimatedStateVector[ 1 ] = -3;

    Eigen::Matrix2d initialEstimatedStateCovarianceMatrix = Eigen::Matrix2d::Zero( );
    initialEstimatedStateCovarianceMatrix( 0, 0 ) = 100;
    initialEstimatedStateCovarianceMatrix( 1, 1 ) = 100;

    // Set system and measurement uncertainty
    Eigen::Matrix2d systemUncertainty = Eigen::Matrix2d::Zero( );
    Eigen::Vector1d measurementUncertainty = Eigen::Vector1d::Zero( );
    systemUncertainty( 0, 0 ) = 100;
    systemUncertainty( 1, 1 ) = 100;
    measurementUncertainty[ 0 ] = 100;

    // Set integrator settings
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            boost::make_shared< numerical_integrators::IntegratorSettings< > > (
                numerical_integrators::euler, initialTime, timeStep );

    // Create control class
    boost::shared_ptr< ControlWrapper< double, double, 2 > > control =
            boost::make_shared< ControlWrapper< double, double, 2 > >( boost::lambda::constant( Eigen::Vector2d::Zero( ) ) );

    // Create extended Kalman filter object
    KalmanFilterDoublePointer extendedFilter = boost::make_shared< ExtendedKalmanFilterDouble >(
                boost::bind( &stateFunction1, _1, _2,
                             boost::bind( &ControlWrapper< double, double, 2 >::getControlVector, control ) ),
                boost::bind( &measurementFunction1, _1, _2 ),
                boost::bind( &stateJacobianFunction1, _1, _2,
                             boost::bind( &ControlWrapper< double, double, 2 >::getControlVector, control ) ),
                boost::lambda::constant( Eigen::Matrix2d::Identity( ) ),
                boost::bind( &measurementJacobianFunction1, _1, _2 ),
                boost::lambda::constant( Eigen::Vector1d::Identity( ) ),
                systemUncertainty, measurementUncertainty,
                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
                integratorSettings );

    // Loop over each time step
    const bool showProgress = false;
    double currentTime = initialTime;
    Eigen::Vector2d currentActualStateVector = initialStateVector;
    Eigen::Vector2d currentNoisyStateVector;
    Eigen::Vector2d currentControlVector = Eigen::Vector2d::Zero( );
    Eigen::Vector1d currentMeasurementVector;
    std::map< double, Eigen::Vector2d > actualStateVectorHistory;
    std::map< double, Eigen::Vector1d > measurementVectorHistory;
    actualStateVectorHistory[ initialTime ] = initialStateVector;
    for( unsigned int i = 0; i < numberOfTimeSteps; i++ )
    {
        // Compute actual values and perturb them
        currentTime += timeStep;
        currentActualStateVector += stateFunction1( currentTime, currentActualStateVector, currentControlVector ) * timeStep;
        currentNoisyStateVector = currentActualStateVector + extendedFilter->produceSystemNoise( ) * timeStep;
        currentMeasurementVector = measurementFunction1( currentTime, currentActualStateVector ) +
                extendedFilter->produceMeasurementNoise( );
        actualStateVectorHistory[ currentTime ] = currentActualStateVector;
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Update control class
        control->setControlVector( currentTime, extendedFilter->getCurrentStateEstimate( ) );

        // Update filter
        extendedFilter->updateFilter( currentTime, currentMeasurementVector );

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "Estimated State: " << extendedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl << std::endl;
        }
    }

    // Check that final state is as expected
    Eigen::Vector2d expectedFinalState = Eigen::Vector2d::Zero( );
    expectedFinalState << 4.5372985336740985, -9.5556424725983025;
    for ( int i = 0; i < expectedFinalState.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( extendedFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ],
                           std::numeric_limits< double >::epsilon( ) );
    }

    // Check that noise is actually normally distributed (within 5 %)
    std::pair< std::vector< Eigen::VectorXd >, std::vector< Eigen::VectorXd > > noiseHistory = extendedFilter->getNoiseHistory( );
    Eigen::MatrixXd systemNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.first );
    Eigen::MatrixXd measurementNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.second );
    for ( unsigned int i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( systemNoise.row( i ) ),
                                    std::sqrt( systemUncertainty( i, i ) ), 5e-2 );
    }
    BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( measurementNoise.row( 0 ) ),
                                std::sqrt( measurementUncertainty( 0, 0 ) ), 5e-2 );

    // Save actual state history
    input_output::writeDataMapToTextFile( actualStateVectorHistory, "EKFActualStateHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF" );

    // Save estimated state history
    input_output::writeDataMapToTextFile( extendedFilter->getEstimatedStateHistory( ), "EKFEstimatedStateHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF" );

    // Save estimated covariance matrix history
    input_output::writeDataMapToTextFile( extendedFilter->getEstimatedCovarianceHistory( ), "EKFEstimatedCovarianceHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF" );

    // Save measurement history
    input_output::writeDataMapToTextFile( measurementVectorHistory, "EKFMeasurementHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KF" );
}

// Constant parameters for example
const double gravitationalParameter = 32.2;

// Functions for extended Kalman filter.
Eigen::Vector3d stateFunction2( const double time, const Eigen::Vector3d& state, const Eigen::Vector3d& control )
{
    Eigen::Vector3d stateDerivative = Eigen::Vector3d::Zero( );
    stateDerivative[ 0 ] = state[ 1 ];
    stateDerivative[ 1 ] = 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            std::pow( state[ 1 ], 2 ) / ( 2.0 * state[ 2 ] ) - gravitationalParameter;
    return stateDerivative;
}
Eigen::Vector1d measurementFunction2( const double time, const Eigen::Vector3d& state )
{
    Eigen::Vector1d measurement;
    measurement[ 0 ] = state[ 0 ];
    return measurement;
}
Eigen::Matrix3d stateJacobianFunction2( const double time, const Eigen::Vector3d& state, const Eigen::Vector3d& control )
{
    Eigen::Matrix3d stateJacobian = Eigen::Matrix3d::Zero( );
    stateJacobian( 0, 1 ) = 1.0;
    stateJacobian( 1, 0 ) = - 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            std::pow( state[ 1 ], 2 ) / ( 44000.0 * state[ 2 ] );
    stateJacobian( 1, 1 ) = 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            state[ 1 ] / state[ 2 ];
    stateJacobian( 1, 2 ) = - 0.0034 * gravitationalParameter * std::exp( - state[ 0 ] / 22000.0 ) *
            std::pow( state[ 1 ], 2 ) / ( 2.0 * std::pow( state[ 2 ], 2 ) );
    return stateJacobian;
}
Eigen::RowVector3d measurementJacobianFunction2( const double time, const Eigen::Vector3d& state )
{
    Eigen::RowVector3d measurementJacobian = Eigen::RowVector3d::Zero( );
    measurementJacobian[ 0 ] = 1.0;
    return measurementJacobian;
}

// Test implementation of extended Kalman filter class.
BOOST_AUTO_TEST_CASE( testExtendedKalmanFilterSecondCase )
{
    using namespace tudat::filters;

    // Set initial conditions
    const double initialTime = 0;
    const double timeStep = 0.1;
    const unsigned int numberOfTimeSteps = 300;

    Eigen::Vector3d initialStateVector;
    initialStateVector[ 0 ] = 200000.0;
    initialStateVector[ 1 ] = -6000.0;
    initialStateVector[ 2 ] = 500.0;

    Eigen::Vector3d initialEstimatedStateVector;
    initialEstimatedStateVector[ 0 ] = 200025.0;
    initialEstimatedStateVector[ 1 ] = -6150.0;
    initialEstimatedStateVector[ 2 ] = 800.0;

    Eigen::Matrix3d initialEstimatedStateCovarianceMatrix = Eigen::Matrix3d::Zero( );
    initialEstimatedStateCovarianceMatrix( 0, 0 ) = std::pow( 1000.0, 2 );
    initialEstimatedStateCovarianceMatrix( 1, 1 ) = 20000.0;
    initialEstimatedStateCovarianceMatrix( 2, 2 ) = std::pow( 300.0, 2 );

    // Set system and measurement uncertainty
    Eigen::Matrix3d systemUncertainty = Eigen::Matrix3d::Zero( );
    Eigen::Vector1d measurementUncertainty = Eigen::Vector1d::Zero( );
    systemUncertainty( 0, 0 ) = std::pow( 100.0, 2 );
    systemUncertainty( 1, 1 ) = std::pow( 10.0, 2 );
    systemUncertainty( 2, 2 ) = std::pow( 1.0, 2 );
    measurementUncertainty[ 0 ] = std::pow( 25.0, 2 );

    // Set integrator settings
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            boost::make_shared< numerical_integrators::IntegratorSettings< > > (
                numerical_integrators::euler, initialTime, timeStep );

    // Create control class
    boost::shared_ptr< ControlWrapper< double, double, 3 > > control =
            boost::make_shared< ControlWrapper< double, double, 3 > >( boost::lambda::constant( Eigen::Vector3d::Zero( ) ) );

    // Create extended Kalman filter object
    KalmanFilterDoublePointer extendedFilter = boost::make_shared< ExtendedKalmanFilterDouble >(
                boost::bind( &stateFunction2, _1, _2,
                             boost::bind( &ControlWrapper< double, double, 3 >::getControlVector, control ) ),
                boost::bind( &measurementFunction2, _1, _2 ),
                boost::bind( &stateJacobianFunction2, _1, _2,
                             boost::bind( &ControlWrapper< double, double, 3 >::getControlVector, control ) ),
                boost::lambda::constant( Eigen::Matrix3d::Zero( ) ),
                boost::bind( &measurementJacobianFunction2, _1, _2 ),
                boost::lambda::constant( Eigen::Vector1d::Zero( ) ),
                systemUncertainty, measurementUncertainty,
                initialTime, initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix,
                integratorSettings );

    // Loop over each time step
    const bool showProgress = false;
    double currentTime = initialTime;
    Eigen::Vector3d currentActualStateVector = initialStateVector;
    Eigen::Vector3d currentNoisyStateVector;
    Eigen::Vector3d currentControlVector = Eigen::Vector3d::Zero( );
    Eigen::Vector1d currentMeasurementVector;
    std::map< double, Eigen::Vector3d > actualStateVectorHistory;
    std::map< double, Eigen::Vector1d > measurementVectorHistory;
    actualStateVectorHistory[ initialTime ] = initialStateVector;
    for( unsigned int i = 0; i < numberOfTimeSteps; i++ )
    {
        // Compute actual values and perturb them
        currentTime += timeStep;
        currentActualStateVector += stateFunction2( currentTime, currentActualStateVector, currentControlVector ) * timeStep;
        currentNoisyStateVector = currentActualStateVector + extendedFilter->produceSystemNoise( ) * timeStep;
        currentMeasurementVector = measurementFunction2( currentTime, currentActualStateVector ) +
                extendedFilter->produceMeasurementNoise( );
        actualStateVectorHistory[ currentTime ] = currentActualStateVector;
        measurementVectorHistory[ currentTime ] = currentMeasurementVector;

        // Update control class
        control->setControlVector( currentTime, extendedFilter->getCurrentStateEstimate( ) );

        // Update filter
        extendedFilter->updateFilter( currentTime, currentMeasurementVector );

        // Print progress
        if ( showProgress )
        {
            std::cout << "Time: " << currentTime << std::endl
                      << "Measurement: " << currentMeasurementVector.transpose( ) << std::endl
                      << "Estimated State: " << extendedFilter->getCurrentStateEstimate( ).transpose( ) << std::endl << std::endl;
        }
    }

    // Check that final state is as expected
    Eigen::Vector3d expectedFinalState = Eigen::Vector3d::Zero( );
    expectedFinalState << 25189.963892316107, -3313.5632678177799, 496.56194076280673;
    for ( int i = 0; i < expectedFinalState.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( extendedFilter->getCurrentStateEstimate( )[ i ] - expectedFinalState[ i ],
                           std::numeric_limits< double >::epsilon( ) );
    }

    // Check that noise is actually normally distributed (within 5 %)
    std::pair< std::vector< Eigen::VectorXd >, std::vector< Eigen::VectorXd > > noiseHistory = extendedFilter->getNoiseHistory( );
    Eigen::MatrixXd systemNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.first );
    Eigen::MatrixXd measurementNoise = utilities::convertStlVectorToEigenMatrix( noiseHistory.second );
    for ( unsigned int i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( systemNoise.row( i ) ),
                                    std::sqrt( systemUncertainty( i, i ) ), 5e-2 );
    }
    BOOST_CHECK_CLOSE_FRACTION( statistics::computeStandardDeviationOfVectorComponents( measurementNoise.row( 0 ) ),
                                std::sqrt( measurementUncertainty( 0, 0 ) ), 5e-2 );

    // Save actual state history
    input_output::writeDataMapToTextFile( actualStateVectorHistory, "EKFActualStateHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook" );

    // Save estimated state history
    input_output::writeDataMapToTextFile( extendedFilter->getEstimatedStateHistory( ), "EKFEstimatedStateHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook" );

    // Save estimated covariance matrix history
    input_output::writeDataMapToTextFile( extendedFilter->getEstimatedCovarianceHistory( ), "EKFEstimatedCovarianceHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook" );

    // Save measurement history
    input_output::writeDataMapToTextFile( measurementVectorHistory, "EKFMeasurementHistory.dat",
                                          "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook" );

    // Save noise histories
    systemNoise.transposeInPlace( );
    measurementNoise.transposeInPlace( );
    input_output::writeMatrixToFile( systemNoise, "systemNoise.dat", 16,
                                     "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook" );
    input_output::writeMatrixToFile( measurementNoise, "measurementNoise.dat", 16,
                                     "/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Test/SimulationOutput/KFBook" );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat