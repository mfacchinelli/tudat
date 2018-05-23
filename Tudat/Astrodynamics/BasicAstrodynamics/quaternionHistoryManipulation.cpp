/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Schaub, H. and Junkins, J., Analytical Mechanics of Space Systems, 2nd ed., ser. Education Series.
 *          American Institute of Aeronautics and Astronautics, 2002.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/quaternionHistoryManipulation.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

#include "Tudat/Mathematics/Statistics/basicStatistics.h"

#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Transform quaternion to opposite rotation.
void convertQuaternionHistoryToMatchSigns( std::map< double, Eigen::Vector4d >& quaternionHistoryMap )
{
    // Get total number of rows
    unsigned int numberOfRows = quaternionHistoryMap.size( );

    // Create concatenation of vectors
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > quaternionHistoryMapSplit =
            utilities::extractKeyAndValuesFromMap< double, double, 4 >( quaternionHistoryMap );
    Eigen::VectorXd timeHistoryVector = quaternionHistoryMapSplit.first;
    Eigen::MatrixXd quaternionHistoryMatrix = quaternionHistoryMapSplit.second;

    // Loop over each quaternion element
    std::vector< unsigned int > indices;
    for ( unsigned int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
    {
        // Compute numerical derivative of quaternion element
        Eigen::VectorXd currentQuaternionElement = quaternionHistoryMatrix.row( quaternionIndex ).transpose( );
        Eigen::VectorXd timeDerivativeOfCurrentQuaternionElement;
        timeDerivativeOfCurrentQuaternionElement.resize( numberOfRows - 1, 1 );
        for ( unsigned int i = 0; i < numberOfRows - 1; i++ )
        {
            timeDerivativeOfCurrentQuaternionElement[ i ] =
                    ( currentQuaternionElement[ i + 1 ] - currentQuaternionElement[ i ] ) /
                    ( timeHistoryVector[ i + 1 ] - timeHistoryVector[ i ] );
        }

        // Compute difference in time derivative
        Eigen::VectorXd differenceInTimeDerivative;
        differenceInTimeDerivative.resize( numberOfRows - 2, 1 );
        for ( unsigned int i = 0; i < numberOfRows - 2; i++ )
        {
            differenceInTimeDerivative[ i ] = timeDerivativeOfCurrentQuaternionElement[ i + 1 ] -
                    timeDerivativeOfCurrentQuaternionElement[ i ];
        }

        // Get standard deviation of first elements and use it as threshold
        double threshold = statistics::computeStandardDeviationOfVectorComponents( differenceInTimeDerivative );

        // Get indices of derivatives beyond the threshold
        for ( unsigned int i = 0; i < differenceInTimeDerivative.size( ); i++ )
        {
            if ( differenceInTimeDerivative[ i ] > threshold )
            {
                indices.push_back( i );
            }
        }
    }

    // Sort indices (duplicates are not considered)
    std::sort( indices.begin( ), indices.end( ) );

    // Invert quaternions based on computed indices
    if ( !indices.empty( ) )
    {
        // Loop over indices
        for ( unsigned int i = 0; i < indices.size( ); i++ )
        {
            if ( i != indices.size( ) - 1 )
            {
                // Skip consecutive indices
                if ( ( indices.at( i + 1 ) - indices.at( i ) ) > 3 )
                {
                    // Invert quaternion
                    for ( unsigned int j = indices.at( i ); j < numberOfRows - 1; j++ )
                    {
                        quaternionHistoryMap[ timeHistoryVector[ j + 1 ] ] *= - 1.0;
                    }
                }
            }
            else
            {
                // Invert quaternion
                for ( unsigned int j = indices.at( i ); j < numberOfRows - 1; j++ )
                {
                    quaternionHistoryMap[ timeHistoryVector[ j + 1 ] ] *= - 1.0;
                }
            }
        }
    }
}

//! Transform quaternion in translational or rotational state to opposite rotation.
void convertQuaternionHistoryToMatchSigns( std::map< double, Eigen::VectorXd >& stateHistoryMap,
                                           const propagators::IntegratedStateType stateType )
{
    // Select index based on input state type
    int quaternionStartIndex;
    switch ( stateType )
    {
    case propagators::translational_state:
        quaternionStartIndex = static_cast< int >( etaUSM7Index );
        break;
    case propagators::rotational_state:
        quaternionStartIndex = static_cast< int >( etaQuaternionIndex );
        break;
    default:
        throw std::runtime_error( "Error in conversion of quaternion history."
                                  "Only translational and rotational propagators are supported." );
    }

    // Extract exponential map from state and convert history
    std::map< double, Eigen::Vector4d > quaternionHistoryMap;
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = stateHistoryMap.begin( );
          mapIterator != stateHistoryMap.end( ); mapIterator++ )
    {
        quaternionHistoryMap[ mapIterator->first ] = mapIterator->second.segment( quaternionStartIndex, 4 );
    }
    convertQuaternionHistoryToMatchSigns( quaternionHistoryMap );

    // Replace state elements with new quaterions
    for ( std::map< double, Eigen::VectorXd >::iterator mapIterator = stateHistoryMap.begin( );
          mapIterator != stateHistoryMap.end( ); mapIterator++ )
    {
        mapIterator->second.segment( quaternionStartIndex, 4 ) = quaternionHistoryMap[ mapIterator->first ];
    }
}

} // namespace orbital_element_conversions

} // namespace tudat