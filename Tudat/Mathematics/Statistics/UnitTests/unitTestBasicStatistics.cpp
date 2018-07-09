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

#include <limits>
#include <iostream>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/Statistics/basicStatistics.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_basic_statistics )

//! Test if sample mean is computed correctly.
BOOST_AUTO_TEST_CASE( testSampleMean )
{
    // Test computation of sample mean on finite population using unbiased estimators.
    // The expected values are computed using the Microsoft Excel the AVERAGE( ) function.

    // Declare vector of sample data.
    std::vector< double > sampleData;

    // Populate vector with sample data.
    sampleData.push_back( 2.5 );
    sampleData.push_back( 6.4 );
    sampleData.push_back( 8.9 );
    sampleData.push_back( 12.7 );
    sampleData.push_back( 15.0 );

    // Set expected sample mean.
    double expectedSampleMean = 9.1;

    // Compute sample mean.
    double computedSampleMean = statistics::computeSampleMean( sampleData );

    // Check if computed sample mean matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSampleMean, expectedSampleMean,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if sample variance is computed correctly.
BOOST_AUTO_TEST_CASE( testSampleVariance )
{
    // Test computation of sample variance on finite population using unbiased estimators.
    // The expected values are computed using the Microsoft Excel the VAR( ) function.

    // Declare vector of sample data.
    std::vector< double > sampleData;

    // Populate vector with sample data.
    sampleData.push_back( 2.5 );
    sampleData.push_back( 6.4 );
    sampleData.push_back( 8.9 );
    sampleData.push_back( 12.7 );
    sampleData.push_back( 15.0 );

    // Declare expected sample variance.
    double expectedSampleVariance = 24.665;

    // Compute sample variance.
    double computedSampleVariance = statistics::computeSampleVariance( sampleData );

    // Check if computed sample variance matches expected value.
    BOOST_CHECK_CLOSE_FRACTION( computedSampleVariance, expectedSampleVariance,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test if moving average is computed correctly. Results compared with MATLAB smooth function.
BOOST_AUTO_TEST_CASE( testMovingAverage )
{
    Eigen::VectorXd vectorOfPoints;
    vectorOfPoints.resize( 15 );
    vectorOfPoints << 0.93329860101525, 3.93372816267124, 1.35032100135611, 2.97099423629127, 3.18245216750598, 1.43494398584927,
            0.915460520182276, 2.60394635060288, 2.09834777464011, 3.04137361348961, 0.265830887303261, 2.96918626998768,
            1.23234701262448, 2.42638755740894, 1.6271912582765;

    Eigen::VectorXd expectedMovingAverage;
    expectedMovingAverage.resize( 15 );
    expectedMovingAverage << 0.93329860101525, 2.07244925501420, 2.47415883376797, 2.57448791073478, 1.97083438223698, 2.22155945208634,
            2.04703015975610, 2.01881444895283, 1.78499182924363, 2.19573697920471, 1.92141711160903, 1.98702506816280, 1.70418859712017,
            1.76197527610331, 1.62719125827650;

    Eigen::VectorXd computedMovingAverage = statistics::computeMovingAverage( vectorOfPoints, 5 );

    // Compare results
    double tolerance = std::numeric_limits< double >::epsilon( );
    for ( unsigned int i = 0; i < 15; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( computedMovingAverage[ i ], expectedMovingAverage[ i ],
                                    tolerance );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
