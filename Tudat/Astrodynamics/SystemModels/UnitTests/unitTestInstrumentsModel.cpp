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

#include "Tudat/Mathematics/Statistics/basicStatistics.h"

using namespace tudat::system_models;

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_instruments_model )

BOOST_AUTO_TEST_CASE( testGeneralModel )
{
    // Declare instrument model
    boost::shared_ptr< InstrumentsModelGenerator > modelGenerator = boost::make_shared< InstrumentsModelGenerator >( );
    boost::shared_ptr< InstrumentsModel > instrumentsModel = modelGenerator->getInstrumentsModel( );

    // Check that nothing works :p
    // Inertial measurement unit
    {
        bool failureDetected = false;
        try
        {
            instrumentsModel->getCurrentInertialMeasurementUnitMeasurement( );
        }
        catch ( const std::runtime_error& caughtException )
        {
            TUDAT_UNUSED_PARAMETER( caughtException );
            failureDetected = true;
        }
        BOOST_CHECK_BITWISE_EQUAL( failureDetected, true );
    }

    // Star tracker
    {
        bool failureDetected = false;
        try
        {
            instrumentsModel->getCurrentStarTrackerMeasurement( );
        }
        catch ( const std::runtime_error& caughtException )
        {
            TUDAT_UNUSED_PARAMETER( caughtException );
            failureDetected = true;
        }
        BOOST_CHECK_BITWISE_EQUAL( failureDetected, true );
    }

    // Deep space network
    {
        bool failureDetected = false;
        try
        {
            instrumentsModel->getCurrentDeepSpaceNetworkMeasurement( );
        }
        catch ( const std::runtime_error& caughtException )
        {
            TUDAT_UNUSED_PARAMETER( caughtException );
            failureDetected = true;
        }
        BOOST_CHECK_BITWISE_EQUAL( failureDetected, true );
    }

    // Altimeter
    {
        bool failureDetected = false;
        try
        {
            instrumentsModel->getCurrentAltimeterMeasurement( );
        }
        catch ( const std::runtime_error& caughtException )
        {
            TUDAT_UNUSED_PARAMETER( caughtException );
            failureDetected = true;
        }
        BOOST_CHECK_BITWISE_EQUAL( failureDetected, true );
    }

    // Generic ranging system
    {
        bool failureDetected = false;
        try
        {
            instrumentsModel->getCurrentGenericRangingSystemMeasurement( );
        }
        catch ( const std::runtime_error& caughtException )
        {
            TUDAT_UNUSED_PARAMETER( caughtException );
            failureDetected = true;
        }
        BOOST_CHECK_BITWISE_EQUAL( failureDetected, true );
    }
}

BOOST_AUTO_TEST_CASE( testInertialMeasurementUnit )
{
    // Declare vectors
    Eigen::Vector3d accelerometerBias, accelerometerScaleFactor, accelerometerNoise, gyroscopeBias, gyroscopeScaleFactor, gyroscopeNoise;
    Eigen::Vector6d accelerometerMisalignment, gyroscopeMisalignment;

    for ( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        // Declare instrument model
        boost::shared_ptr< InstrumentsModelGenerator > modelGenerator = boost::make_shared< InstrumentsModelGenerator >( );
        boost::shared_ptr< InstrumentsModel > instrumentsModel = modelGenerator->getInstrumentsModel( );

        // Set accuracies
        switch ( testCase )
        {
        case 0: // without errors
        {
            // Set accelerometer accuracies
            accelerometerBias.setZero( );
            accelerometerScaleFactor.setZero( );
            accelerometerMisalignment.setZero( );
            accelerometerNoise.setZero( );

            // Set gyroscope accuracies
            gyroscopeBias.setZero( );
            gyroscopeScaleFactor.setZero( );
            gyroscopeMisalignment.setZero( );
            gyroscopeNoise.setZero( );
            break;
        }
        case 1: // with noise
        {
            // Set accelerometer noise
            accelerometerNoise << 0.01, 0.005, 0.2;

            // Set gyroscope noise
            gyroscopeNoise << 0.1, 0.05, 0.002;
            break;
        }
        case 2: // with errors
        {
            // Set accelerometer accuracies
            accelerometerBias = Eigen::Vector3d::Constant( -1.2345e-2 );
            accelerometerScaleFactor = Eigen::Vector3d::Constant( 5.4321e-3 );
            accelerometerMisalignment = Eigen::Vector6d::Constant( 3.4152e-6 );
            accelerometerNoise.setZero( );

            // Set gyroscope accuracies
            gyroscopeBias = Eigen::Vector3d::Constant( 3.4152e-2 );
            gyroscopeScaleFactor = Eigen::Vector3d::Constant( 1.2345e-3 );
            gyroscopeMisalignment = Eigen::Vector6d::Constant( -5.4321e-6 );
            gyroscopeNoise.setZero( );
            break;
        }
        }

        // Add inertial measurement unit
        instrumentsModel->addInertialMeasurementUnit( accelerometerBias, accelerometerScaleFactor, accelerometerMisalignment, accelerometerNoise,
                                                      gyroscopeBias, gyroscopeScaleFactor, gyroscopeMisalignment, gyroscopeNoise );

        // Check that first value is indeed NaN (instrument model has not been updated)
        if ( testCase == 0 )
        {
            for ( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_EQUAL( true, std::isnan( instrumentsModel->getCurrentAccelerometerMeasurement( )[ i ] ) );
                BOOST_CHECK_EQUAL( true, std::isnan( instrumentsModel->getCurrentGyroscopeMeasurement( )[ i ] ) );
            }
        }

        // Update environment
        modelGenerator->updateEnvironment( );

        // Check measured IMU data
        double expectedAccelerationMagnitude = 1.0e-8;
        switch ( testCase )
        {
        case 0:
        {
            // Accelerometer measurement should return solar radiation, which for Mars, is about 1e-8 in magnitude
            BOOST_CHECK_SMALL( std::fabs( expectedAccelerationMagnitude -
                                          instrumentsModel->getCurrentAccelerometerMeasurement( ).norm( ) ), 5.0e-8 );

            // Gyroscope measurment should be zero, since no rotational motion is implemented
            BOOST_CHECK_EQUAL( Eigen::Vector3d::Zero( ), instrumentsModel->getCurrentGyroscopeMeasurement( ) );
            break;
        }
        case 2:
        {
            // Accelerometer measurement should return solar radiation, which for Mars, is about 1e-8 in magnitude
            // However, the value of errors are also to be included
            expectedAccelerationMagnitude += accelerometerBias.norm( );
            BOOST_CHECK_SMALL( expectedAccelerationMagnitude - instrumentsModel->getCurrentAccelerometerMeasurement( ).norm( ), 1.0e-1 );

            // Gyroscope measurment should be zero, since no rotational motion is implemented
            // However, the value of errors are also to be included
            BOOST_CHECK_EQUAL( gyroscopeBias.norm( ), instrumentsModel->getCurrentGyroscopeMeasurement( ).norm( ) );
            break;
        }
        case 1:
        {
            // Propagate for a few minutes and store results, to check that noise is within STD
            std::vector< Eigen::VectorXd > vectorOfGyroscopeMeasurements;
            vectorOfGyroscopeMeasurements.resize( 1000 );
            for ( unsigned int i = 0; i < 1000; i++ )
            {
                // Update environment
                modelGenerator->updateEnvironment( );

                // Store measurement
                vectorOfGyroscopeMeasurements.at( i ) = instrumentsModel->getCurrentGyroscopeMeasurement( );
            }

            // Compute STD
            Eigen::VectorXd varianceOfMeasurements = statistics::computeSampleVariance( vectorOfGyroscopeMeasurements );

            // Verify STD
            for ( unsigned i = 0; i < 3; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( gyroscopeNoise[ i ], std::sqrt( varianceOfMeasurements[ i ] ), 5.0e-2 );
            }
        }
        }
    }
}

//BOOST_AUTO_TEST_CASE( testStarTracker )
//{
//    // Declare instrument model
//    boost::shared_ptr< InstrumentsModelGenerator > modelGenerator = boost::make_shared< InstrumentsModelGenerator >( );
//    boost::shared_ptr< InstrumentsModel > instrumentsModel = modelGenerator->getInstrumentsModel( );

//}

//BOOST_AUTO_TEST_CASE( testDeepSpaceNetwork )
//{
//    // Declare instrument model
//    boost::shared_ptr< InstrumentsModelGenerator > modelGenerator = boost::make_shared< InstrumentsModelGenerator >( );
//    boost::shared_ptr< InstrumentsModel > instrumentsModel = modelGenerator->getInstrumentsModel( );

//}

//BOOST_AUTO_TEST_CASE( testAltimeter )
//{
//    // Declare instrument model
//    boost::shared_ptr< InstrumentsModelGenerator > modelGenerator = boost::make_shared< InstrumentsModelGenerator >( );
//    boost::shared_ptr< InstrumentsModel > instrumentsModel = modelGenerator->getInstrumentsModel( );

//}

//BOOST_AUTO_TEST_CASE( testGenericRangingSystem )
//{
//    // Declare instrument model
//    boost::shared_ptr< InstrumentsModelGenerator > modelGenerator = boost::make_shared< InstrumentsModelGenerator >( );
//    boost::shared_ptr< InstrumentsModel > instrumentsModel = modelGenerator->getInstrumentsModel( );

//}

BOOST_AUTO_TEST_SUITE_END( )

}

}

