/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MULTI_DIMENSIONAL_ARRAY_WRITER_H
#define TUDAT_MULTI_DIMENSIONAL_ARRAY_WRITER_H

#include <map>

#include <boost/multi_array.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace input_output
{

//! Function to print to console which aerodynamic coefficients are being saved.
/*!
 *  Function to print to console which aerodynamic coefficients are being saved.
 *  \param index Index of coefficient to be saved.
 */
void informUserOnSavedCoefficient( const int index );

//! Interface class for writing coefficients as a function of N independent variables to a file.
/*!
 *  Interface class for writing coefficients as a function of N independent variables to a file. This class is used instead
 *  of a single templated free function to allow multi-arrays of different sizes to be created using the same interface.
 *  NOTE: The possibility of using a single templated implementation for arbitrary multi-array size should be investigated
 *  in the future.
 */
template< unsigned int NumberOfIndependentVariables, unsigned int NumberOfCoefficients >
class MultiArrayFileWriter
{
public:

    //! Function to write one-dimensional array data and independent variables to file.
    /*!
     *  Function to write one-dimensional array data and independent variables to file.
     *  \tparam NumberOfIndependentVariables Number of independent variables.
     *  \tparam NumberOfCoefficients Number of coefficients.
     *  \param fileName File where coefficients need to be saved. All coefficients are saved in the same file, with the
     *      independent variable in the first column.
     *  \param coefficientIndices Indices of coefficients to be saved to file.
     *  \param independentVariables Vector of independent variable vectors.
     *  \param dependentVariables Multi-array of dependent variables for each independent variable.
     */
    static void writeMultiArrayAndIndependentVariablesToFiles( const std::string& fileName,
                                                               const std::vector< int >& coefficientIndices,
                                                               const std::vector< std::vector< double > >& independentVariables,
                                                               const boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients,
                                                               1 >, NumberOfIndependentVariables >& dependentVariables );

    //! Function to write multi-dimensional array data and independent variables to file.
    /*!
     *  Function to write multi-dimensional array data and independent variables to file.
     *  \tparam NumberOfIndependentVariables Number of independent variables.
     *  \tparam NumberOfCoefficients Number of coefficients.
     *  \param fileNamesMap Map of files where coefficients need to be saved. Each coefficient is saved in a separate file, thus
     *      number of files has to match number of coefficients.
     *  \param independentVariables Vector of independent variable vectors.
     *  \param dependentVariables Multi-array of dependent variables for each independent variable.
     */
    static void writeMultiArrayAndIndependentVariablesToFiles( const std::map< int, std::string >& fileNamesMap,
                                                               const std::vector< std::vector< double > >& independentVariables,
                                                               const boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients,
                                                               1 >, NumberOfIndependentVariables >& dependentVariables );

};

//! Interface class for writing coefficients as a function of 1 independent variables to a file.
template< unsigned int NumberOfCoefficients >
class MultiArrayFileWriter< 1, NumberOfCoefficients >
{
public:

    //! Function to write one-dimensional array data and independent variables to file.
    /*!
     *  Function to write one-dimensional array data and independent variables to file.
     *  \tparam NumberOfIndependentVariables Number of independent variables.
     *  \tparam NumberOfCoefficients Number of coefficients.
     *  \param fileName File where coefficients need to be saved. All coefficients are saved in the same file, with the
     *      independent variable in the first column.
     *  \param coefficientIndices Indices of coefficients to be saved to file.
     *  \param independentVariables Vector of independent variable vectors.
     *  \param dependentVariables Multi-array of dependent variables for each independent variable.
     */
    static void writeMultiArrayAndIndependentVariablesToFiles( const std::string& fileName,
                                                               const std::vector< int >& coefficientIndices,
                                                               const std::vector< std::vector< double > >& independentVariables,
                                                               const boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients,
                                                               1 >, 1 >& dependentVariables )
    {
        // Get actual number of coefficients
        // The user may only want to save some coefficients (e.g., (s)he may want to disregard CS)
        unsigned int actualNumberOfCoefficients = coefficientIndices.size( );

        // Check consistency of inputs
        // Check that number of coefficients is not too large
        if ( actualNumberOfCoefficients > NumberOfCoefficients )
        {
            throw std::runtime_error( "Error while saving one-dimensional aerodynamic coefficients to file. "
                                      "Number of coefficients does not match number of input files." );
        }
        // Check that number of independent variables matches vector size
        if ( independentVariables.size( ) != 1 )
        {
            throw std::runtime_error( "Error while saving one-dimensional aerodynamic coefficients to file. "
                                      "Number of independent variables does not match size of independent variables vector." );
        }

        // Inform user on which variable is being saved
        for ( int index: coefficientIndices )
        {
            informUserOnSavedCoefficient( index );
        }

        // Open file
        FILE * fileIdentifier = std::fopen( fileName.c_str( ), "w" );

        // Loop over independent variable rows
        for ( unsigned int j = 0; j < independentVariables.at( 0 ).size( ); j++ )
        {
            // Save independent variable
            fprintf( fileIdentifier, "%.10f ", independentVariables.at( 0 ).at( j ) );

            // Save dependent variables
            for ( unsigned int i = 0; i < actualNumberOfCoefficients; i++ )
            {
                fprintf( fileIdentifier, "%.10f ", dependentVariables[ j ][ coefficientIndices.at( i ) ] );
            }
            fprintf( fileIdentifier, "\n" );
        }

        // Close file
        std::fclose( fileIdentifier );
    }

    //! Function to write one-dimensional array data and independent variables to file.
    /*!
     *  Function to write one-dimensional array data and independent variables to file.
     *  \tparam NumberOfIndependentVariables Number of independent variables.
     *  \tparam NumberOfCoefficients Number of coefficients.
     *  \param fileNamesMap Map of files where coefficients need to be saved. Each coefficient is saved in a separate file, thus
     *      number of files has to match number of coefficients.
     *  \param independentVariables Vector of independent variable vectors.
     *  \param dependentVariables Multi-array of dependent variables for each independent variable.
     */
    static void writeMultiArrayAndIndependentVariablesToFiles( const std::map< int, std::string >& fileNamesMap,
                                                               const std::vector< std::vector< double > >& independentVariables,
                                                               const boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients,
                                                               1 >, 1 >& dependentVariables )
    {
        // Get actual number of coefficients
        // The user may only want to save some coefficients (e.g., (s)he may want to disregard CS)
        unsigned int actualNumberOfCoefficients = fileNamesMap.size( );

        // Check consistency of inputs
        // Check that number of coefficients is not too large
        if ( actualNumberOfCoefficients > NumberOfCoefficients )
        {
            throw std::runtime_error( "Error while saving multi-dimensional aerodynamic coefficients to file. "
                                      "Number of coefficients does not match number of input files." );
        }
        // Check that number of independent variables matches vector size
        if ( independentVariables.size( ) != 1 )
        {
            throw std::runtime_error( "Error while saving multi-dimensional aerodynamic coefficients to file. "
                                      "Number of independent variables does not match size of independent variables vector." );
        }

        // Loop over each file
        FILE * fileIdentifier;
        for ( std::map< int, std::string >::const_iterator fileIterator = fileNamesMap.begin( );
              fileIterator != fileNamesMap.end( ); fileIterator++ )
        {
            // Inform user on which variable is being saved
            informUserOnSavedCoefficient( fileIterator->first );

            // Open file
            fileIdentifier = std::fopen( fileIterator->second.c_str( ), "w" );

            // Print number of independent variables
            fprintf( fileIdentifier, "%d\n\n", 1 );

            // Print independent variables
            for ( unsigned int j = 0; j < 1; j++ )
            {
                for ( unsigned int k = 0; k < independentVariables.at( j ).size( ); k++ )
                {
                    fprintf( fileIdentifier, "%.10f ", independentVariables.at( j ).at( k ) );
                }
                fprintf( fileIdentifier, "\n" );
            }
            fprintf( fileIdentifier, "\n" );

            // Print dependent variables
            for ( unsigned int k = 0; k < independentVariables.at( 0 ).size( ); k++ )
            {
                fprintf( fileIdentifier, "%.10f ", dependentVariables[ k ][ fileIterator->first ] );
                fprintf( fileIdentifier, "\n" );
            }

            // Close file
            std::fclose( fileIdentifier );
        }
    }

};

//! Interface class for writing coefficients as a function of 2 independent variables to a file.
template< unsigned int NumberOfCoefficients >
class MultiArrayFileWriter< 2, NumberOfCoefficients >
{
public:

    //! Function to write multi-dimensional array data and independent variables to file.
    /*!
     *  Function to write multi-dimensional array data and independent variables to file.
     *  \tparam NumberOfIndependentVariables Number of independent variables.
     *  \tparam NumberOfCoefficients Number of coefficients.
     *  \param fileNamesMap Map of files where coefficients need to be saved. Each coefficient is saved in a separate file, thus
     *      number of files has to match number of coefficients.
     *  \param independentVariables Vector of independent variable vectors.
     *  \param dependentVariables Multi-array of dependent variables for each independent variable.
     */
    static void writeMultiArrayAndIndependentVariablesToFiles( const std::map< int, std::string >& fileNamesMap,
                                                               const std::vector< std::vector< double > >& independentVariables,
                                                               const boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients,
                                                               1 >, 2 >& dependentVariables )
    {
        // Get actual number of coefficients
        // The user may only want to save some coefficients (e.g., (s)he may want to disregard CS)
        unsigned int actualNumberOfCoefficients = fileNamesMap.size( );

        // Check consistency of inputs
        // Check that number of coefficients is not too large
        if ( actualNumberOfCoefficients > NumberOfCoefficients )
        {
            throw std::runtime_error( "Error while saving multi-dimensional aerodynamic coefficients to file. "
                                      "Number of coefficients does not match number of input files." );
        }
        // Check that number of independent variables matches vector size
        if ( independentVariables.size( ) != 2 )
        {
            throw std::runtime_error( "Error while saving multi-dimensional aerodynamic coefficients to file. "
                                      "Number of independent variables does not match size of independent variables vector." );
        }

        // Loop over each file
        FILE * fileIdentifier;
        for ( std::map< int, std::string >::const_iterator fileIterator = fileNamesMap.begin( );
              fileIterator != fileNamesMap.end( ); fileIterator++ )
        {
            // Inform user on which variable is being saved
            informUserOnSavedCoefficient( fileIterator->first );

            // Open file
            fileIdentifier = std::fopen( fileIterator->second.c_str( ), "w" );

            // Print number of independent variables
            fprintf( fileIdentifier, "%d\n\n", 2 );

            // Print independent variables
            for ( unsigned int j = 0; j < 2; j++ )
            {
                for ( unsigned int k = 0; k < independentVariables.at( j ).size( ); k++ )
                {
                    fprintf( fileIdentifier, "%.10f ", independentVariables.at( j ).at( k ) );
                }
                fprintf( fileIdentifier, "\n" );
            }
            fprintf( fileIdentifier, "\n" );

            // Print dependent variables
            for ( unsigned int k = 0; k < independentVariables.at( 0 ).size( ); k++ )
            {
                for ( unsigned int l = 0; l < independentVariables.at( 1 ).size( ); l++ )
                {
                    fprintf( fileIdentifier, "%.10f ", dependentVariables[ k ][ l ][ fileIterator->first ] );
                }
                fprintf( fileIdentifier, "\n" );
            }

            // Close file
            std::fclose( fileIdentifier );
        }
    }

};

//! Interface class for writing coefficients as a function of 3 independent variables to a file.
template< unsigned int NumberOfCoefficients >
class MultiArrayFileWriter< 3, NumberOfCoefficients >
{
public:

    //! Function to write multi-dimensional array data and independent variables to file.
    /*!
     *  Function to write multi-dimensional array data and independent variables to file.
     *  \tparam NumberOfIndependentVariables Number of independent variables.
     *  \tparam NumberOfCoefficients Number of coefficients.
     *  \param fileNamesMap Map of files where coefficients need to be saved. Each coefficient is saved in a separate file, thus
     *      number of files has to match number of coefficients.
     *  \param independentVariables Vector of independent variable vectors.
     *  \param dependentVariables Multi-array of dependent variables for each independent variable.
     */
    static void writeMultiArrayAndIndependentVariablesToFiles( const std::map< int, std::string >& fileNamesMap,
                                                               const std::vector< std::vector< double > >& independentVariables,
                                                               const boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients,
                                                               1 >, 3 >& dependentVariables )
    {
        // Get actual number of coefficients
        // The user may only want to save some coefficients (e.g., (s)he may want to disregard CS)
        unsigned int actualNumberOfCoefficients = fileNamesMap.size( );

        // Check consistency of inputs
        // Check that number of coefficients is not too large
        if ( actualNumberOfCoefficients > NumberOfCoefficients )
        {
            throw std::runtime_error( "Error while saving multi-dimensional aerodynamic coefficients to file. "
                                      "Number of coefficients does not match number of input files." );
        }
        // Check that number of independent variables matches vector size
        if ( independentVariables.size( ) != 3 )
        {
            throw std::runtime_error( "Error while saving multi-dimensional aerodynamic coefficients to file. "
                                      "Number of independent variables does not match size of independent variables vector." );
        }

        // Loop over each file
        FILE * fileIdentifier;
        for ( std::map< int, std::string >::const_iterator fileIterator = fileNamesMap.begin( );
              fileIterator != fileNamesMap.end( ); fileIterator++ )
        {
            // Inform user on which variable is being saved
            informUserOnSavedCoefficient( fileIterator->first );

            // Open file
            fileIdentifier = std::fopen( fileIterator->second.c_str( ), "w" );

            // Print number of independent variables
            fprintf( fileIdentifier, "%d\n\n", 3 );

            // Print independent variables
            for ( unsigned int j = 0; j < 3; j++ )
            {
                for ( unsigned int k = 0; k < independentVariables.at( j ).size( ); k++ )
                {
                    fprintf( fileIdentifier, "%.10f ", independentVariables.at( j ).at( k ) );
                }
                fprintf( fileIdentifier, "\n" );
            }
            fprintf( fileIdentifier, "\n" );

            // Print dependent variables
            for ( unsigned int j = 0; j < independentVariables.at( 2 ).size( ); j++ )
            {
                for ( unsigned int k = 0; k < independentVariables.at( 0 ).size( ); k++ )
                {
                    for ( unsigned int l = 0; l < independentVariables.at( 1 ).size( ); l++ )
                    {
                        fprintf( fileIdentifier, "%.10f ", dependentVariables[ k ][ l ][ j ][ fileIterator->first ] );
                    }
                    fprintf( fileIdentifier, "\n" );
                }
                fprintf( fileIdentifier, "\n" );
            }

            // Close file
            std::fclose( fileIdentifier );
        }
    }

};

} // namespace input_output

} // namespace tudat

#endif // TUDAT_MULTI_DIMENSIONAL_ARRAY_WRITER_H