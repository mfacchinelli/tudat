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

#include "Tudat/External/SpartaInterface/spartaInterface.h"

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientGenerator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

namespace tudat
{

namespace sparta_interface
{

// Inherit independent variables and vehicle characteristics from aerodynamic coefficient generator
using aerodynamics::AerodynamicCoefficientGenerator;

//! Function to sort the rows of a matrix, based on the specified column and specified order.
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > sortMatrixRows(
        const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >& matrixToBeSorted,
        const int referenceColumn, const bool descendingOrder )
{
    // Declare eventual output vector
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > sortedMatrix;
    sortedMatrix.resizeLike( matrixToBeSorted );

    // Loop over rows and assign to new matrix
    Eigen::Matrix< double, 1, Eigen::Dynamic, Eigen::RowMajor > currentRow;
    for ( unsigned int i = 0; i < matrixToBeSorted.rows( ); i++ )
    {
        // Retrieve each row in 'scrambled order' and assign it to the new matrix
        currentRow = matrixToBeSorted.row( i );

        // Based on requested order
        if ( descendingOrder )
        {
            sortedMatrix.row( matrixToBeSorted.rows( ) - currentRow[ referenceColumn ] ) = currentRow;
            // no (- 1)'s are needed since both are defined starting from 1
        }
        else
        {
            sortedMatrix.row( currentRow[ referenceColumn ] - 1 ) = currentRow;
        }
    }

    // Give output
    return sortedMatrix;
}

//! Open and read geometry file for SPARTA simulation.
void SpartaInterface::analyzeGeometryFile( const std::string& geometryFileUser,
                                           const double referenceArea,
                                           const Eigen::Vector3d& momentReferencePoint )
{
    // Copy input shape file to default name
    std::string commandString = "cp " + geometryFileUser + " " + input_output::getSpartaInternalGeometryFile( );
    std::system( commandString.c_str( ) );

    // Extract information on vehicle geometry
    std::pair< Eigen::Matrix< double, Eigen::Dynamic, 3 >, Eigen::Matrix< int, Eigen::Dynamic, 3 > >
            geometryData = input_output::readSpartaGeometryFile( geometryFileUser );
    shapePoints_ = geometryData.first;
    shapeTriangles_ = geometryData.second;
    numberOfPoints_ = shapePoints_.rows( );
    numberOfTriangles_ = shapeTriangles_.rows( );

    // Resize pressure and shear matrices
    meanPressureValues_.resize( 3, numberOfTriangles_ );
    meanShearValues_.resize( 3, numberOfTriangles_ );

    // Get maximum and minimum values in each dimension
    maximumDimensions_ = shapePoints_.colwise( ).maxCoeff( );
    minimumDimensions_ = shapePoints_.colwise( ).minCoeff( );

    // Compute normal to surface elements, area of surface elements and moment arm values
    Eigen::Matrix3d currentVertices;
    Eigen::Vector3d currentNormal;
    Eigen::Vector3d currentCentroid;
    double currentNormalNorm;
    elementSurfaceNormal_.resize( 3, numberOfTriangles_ );
    elementSurfaceArea_.resize( 1, numberOfTriangles_ );
    elementMomentArm_.resize( 3, numberOfTriangles_ );
    for ( int i = 0; i < numberOfTriangles_; i++ )
    {
        // Compute properties of current surface element
        for ( unsigned int j = 0; j < 3; j++ )
        {
            currentVertices.row( j ) = shapePoints_.row( shapeTriangles_( i, j ) - 1 );
        }
        currentNormal = ( currentVertices.row( 1 ) - currentVertices.row( 0 ) ).cross(
                    currentVertices.row( 2 ) - currentVertices.row( 0 ) );
        currentNormalNorm = currentNormal.norm( );
        currentCentroid = currentVertices.colwise( ).sum( ) / 3.0;

        // Find normal, area and distance to reference point
        elementSurfaceNormal_.col( i ) = currentNormal / currentNormalNorm;
        elementSurfaceArea_( i ) = 0.5 * currentNormalNorm;
        elementMomentArm_.col( i ) = currentCentroid - momentReferencePoint;
    }

    // Compute cross-sectional area
    for ( unsigned int i = 0; i < 3; i++ )
    {
        shapeCrossSectionalArea_( i ) =
                0.5 * elementSurfaceNormal_.row( i ).cwiseAbs( ).dot( elementSurfaceArea_ );
    }

    // Check consistency with input dimensions
    const double tolerance = 1e-5;
    if ( std::fabs( shapeCrossSectionalArea_( referenceDimension_ ) - referenceArea ) > tolerance )
    {
        throw std::runtime_error( "Error in SPARTA geometry file. Input reference area does not match the "
                                  "combination of reference axis and geometry. Note that the first dimension "
                                  "is identified with 0. Tolerance set to: " + std::to_string( tolerance ) );
    }
}

//! Retrieve simulation conditions based on input and geometry.
void SpartaInterface::getSimulationConditions( )
{
    // Simulation boundary and grid
    for ( unsigned int i = 0; i < 3; i++ )
    {
        simulationBoundaries_( 2 * i ) = minimumDimensions_( i ) + 0.5 * minimumDimensions_.minCoeff( ); // add extra space around shape
        simulationBoundaries_( 2 * i + 1 ) = maximumDimensions_( i ) + 0.5 * maximumDimensions_.maxCoeff( ); // add extra space around shape
        if ( i == referenceDimension_ )
        {
            simulationBoundaries_( 2 * i ) -= 1.0; // add extra space along axis of velocity
            simulationBoundaries_( 2 * i + 1 ) += 1.0; // add extra space along axis of velocity
        }
        simulationGrid_( i ) = simulationBoundaries_( 2 * i + 1 ) - simulationBoundaries_( 2 * i );
    }
    simulationGrid_ /= gridSpacing_;

    // Convert molecular speed ratio to stream velocity and compute simulation time step and ratio of real to simulated variables
    freeStreamVelocities_.resize( dataPointsOfIndependentVariables_.at( 0 ).size( ), dataPointsOfIndependentVariables_.at( 1 ).size( ) );
    simulationTimeStep_.resize( dataPointsOfIndependentVariables_.at( 0 ).size( ), dataPointsOfIndependentVariables_.at( 1 ).size( ) );
    ratioOfRealToSimulatedParticles_.resize( dataPointsOfIndependentVariables_.at( 0 ).size( ), 1 );
    double simulationBoxLengthAlongReferenceAxis = ( simulationBoundaries_( 2 * referenceDimension_ + 1 ) -
                                                     simulationBoundaries_( 2 * referenceDimension_ ) );
    for ( unsigned int h = 0; h < dataPointsOfIndependentVariables_.at( 0 ).size( ); h++ )
    {
        for ( unsigned int m = 0; m < dataPointsOfIndependentVariables_.at( 1 ).size( ); m++ )
        {
            freeStreamVelocities_( h, m ) = dataPointsOfIndependentVariables_.at( 1 ).at( m ) *
                    atmosphericConditions_[ speed_of_sound_index ].at( h );
            simulationTimeStep_( h, m ) = 0.1 * simulationBoxLengthAlongReferenceAxis /
                    freeStreamVelocities_( h, m );
            // time step is taken as time it takes for a particle to travel for 10 % of the box
        }
        ratioOfRealToSimulatedParticles_( h ) = atmosphericConditions_[ number_density_index ].at( h ) *
                std::pow( gridSpacing_, 3 ) / simulatedParticlesPerCell_;
    }
}

//! Run SPARTA simulation.
void SpartaInterface::runSpartaSimulation( const unsigned int h, const unsigned int m, const unsigned int a )
{
    // Generate command string for SPARTA
    std::string runSPARTACommandString = "cd " + input_output::getSpartaDataPath( ) + "; ";
    if ( MPIExecutable_ != "" )
    {
        if ( numberOfCores_ < 1 )
        {
            throw std::runtime_error( "Error in SPARTA rarefied flow analysis. Number of cores needs to be "
                                      "an integer value larger or equal to one." );
        }
        runSPARTACommandString = runSPARTACommandString + MPIExecutable_ +
                " -np " + std::to_string( numberOfCores_ ) + " ";
    }
    runSPARTACommandString = runSPARTACommandString + SPARTAExecutable_ + " -echo log ";
    if ( !printProgressInCommandWindow_ )
    {
        runSPARTACommandString = runSPARTACommandString + "-screen none ";
    }
    runSPARTACommandString = runSPARTACommandString + "-in " + input_output::getSpartaInputFile( );

    // Get velocity vector
    velocityVector_ = Eigen::Vector3d::Zero( );
    velocityVector_( referenceDimension_ ) = ( ( referenceAxis_ >= 0 ) ? - 1.0 : 1.0 ) *
            freeStreamVelocities_( h, m );

    // Print to file (each line corresponds to a different line in the file)
    FILE * fileIdentifier = std::fopen( input_output::getSpartaInputFile( ).c_str( ), "w" );
    std::fprintf( fileIdentifier, inputTemplate_.c_str( ), // user-defined variables below
                  simulationBoundaries_( 0 ), simulationBoundaries_( 1 ), simulationBoundaries_( 2 ), // continues on next line
                  simulationBoundaries_( 3 ), simulationBoundaries_( 4 ), simulationBoundaries_( 5 ),
                  simulationGrid_( 0 ), simulationGrid_( 1 ), simulationGrid_( 2 ),
                  atmosphericConditions_[ number_density_index ].at( h ), ratioOfRealToSimulatedParticles_( h ),
                  simulationGases_.c_str( ),
                  simulationGases_.c_str( ), velocityVector_( 0 ), velocityVector_( 1 ), velocityVector_( 2 ),
                  simulationGases_.c_str( ), atmosphericConditions_[ temperature_index ].at( h ),
                  unit_conversions::convertRadiansToDegrees( dataPointsOfIndependentVariables_.at( 2 ).at( a ) ),
                  wallTemperature_, accommodationCoefficient_,
                  simulationTimeStep_( h, m ),
                  input_output::getSpartaOutputDirectory( ).c_str( ) );
    std::fclose( fileIdentifier );

    // Run SPARTA
    systemStatus_ = std::system( runSPARTACommandString.c_str( ) );
    if ( systemStatus_ != 0 )
    {
        throw std::runtime_error( "Error in SPARTA rarefied flow analysis. "
                                  "SPARTA Simulation failed. See the log.sparta file in "
                                  "Tudat/External/SPARTA/ for more details." );
    }
}

//! Process SPARTA output
Eigen::Vector6d SpartaInterface::processSpartaOutput( const double referenceArea, const double referenceLength,
                                                      const unsigned int h, const unsigned int m )
{
    // Read output files and compute mean pressure and shear force values
    meanPressureValues_.setZero( );
    meanShearValues_.setZero( );
    for ( unsigned int i = 0; i < outputFileExtensions_.size( ); i++ )
    {
        // Read file
        outputMatrix_ = input_output::readMatrixFromFile( temporaryOutputFile_ + outputFileExtensions_.at( i ),
                                                          "\t ;,", "%", 9 );

        // Sort rows, since with multiple cores running, row order is scrambled
        outputMatrix_ = sortMatrixRows( outputMatrix_, 0 );

        // Store pressure and shear values
        for ( unsigned int j = 0; j < 3; j++ )
        {
            meanPressureValues_.row( j ) += outputMatrix_.col( j + 1 ).transpose( );
            meanShearValues_.row( j ) += outputMatrix_.col( j + 4 ).transpose( );
        }
    }

    // Compute mean
    meanPressureValues_ /= outputFileExtensions_.size( );
    meanShearValues_ /= outputFileExtensions_.size( );

    // Clean up results folder
    std::string commandString = "rm " + input_output::getSpartaOutputPath( ) + "/coeff.*";
    std::system( commandString.c_str( ) );

    // Convert pressure and shear forces to aerodynamic coefficients
    return aerodynamics::computeAerodynamicCoefficientsFromPressureShearForces(
                meanPressureValues_,
                meanShearValues_,
                atmosphericConditions_[ density_index ].at( h ),
                atmosphericConditions_[ pressure_index ].at( h ),
                freeStreamVelocities_( h, m ),
                elementSurfaceNormal_,
                elementSurfaceArea_,
                elementMomentArm_,
                referenceArea,
                referenceLength );
}

} // namespace sparta_interface

} // namespace tudat
