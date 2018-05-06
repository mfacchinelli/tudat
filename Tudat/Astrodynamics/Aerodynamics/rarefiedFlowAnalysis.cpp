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
 *      Klothakis, A. and Nikolos, I., “Modeling of Rarefied Hypersonic Flows Using the Massively
 *        Parallel DSMC Kernel “SPARTA”,” in 8th GRACM International Congress on Computational Mechanics,
 *        Volos, Greece, July 2015.
 *      Dirkx, D. and Mooij, E., Conceptual Shape Optimization of Entry Vehicles. Springer, 2017.
 *
 */

#include <boost/make_shared.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Geometry>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/InputOutput/spartaDataReader.h"
#include "Tudat/InputOutput/spartaInputOutput.h"

#include "Tudat/Astrodynamics/Aerodynamics/rarefiedFlowAnalysis.h"

namespace tudat
{

namespace aerodynamics
{

using namespace unit_conversions;
using namespace input_output;

//! Returns default values of number density for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowNumberDensityPoints( )
{
    std::vector< double > numberDensityPoints;
    numberDensityPoints.resize( 5 );
    numberDensityPoints[ 0 ] = 125.0e3;
    numberDensityPoints[ 1 ] = 150.0e3;
    numberDensityPoints[ 2 ] = 200.0e3;
    numberDensityPoints[ 3 ] = 300.0e3;
    numberDensityPoints[ 4 ] = 500.0e3;

    return numberDensityPoints;
}

//! Returns default values of molecular speed ratio for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowMolecularSpeedRatioPoints(
        const std::string& molecularSpeedRatioRegime )
{
    std::vector< double > molecularSpeedRatioPoints;

    // Set default points for full speed analysis.
    if ( molecularSpeedRatioRegime == "Full" )
    {
        molecularSpeedRatioPoints.resize( 7 );
        molecularSpeedRatioPoints[ 0 ] = 1.0;
        molecularSpeedRatioPoints[ 1 ] = 2.5;
        molecularSpeedRatioPoints[ 2 ] = 5.0;
        molecularSpeedRatioPoints[ 3 ] = 10.0;
        molecularSpeedRatioPoints[ 4 ] = 25.0;
        molecularSpeedRatioPoints[ 5 ] = 50.0;
        molecularSpeedRatioPoints[ 6 ] = 100.0;
    }
    // Set default points for low speed analysis.
    else if ( molecularSpeedRatioRegime == "Low" )
    {
        molecularSpeedRatioPoints.resize( 4 );
        molecularSpeedRatioPoints[ 0 ] = 1.0;
        molecularSpeedRatioPoints[ 1 ] = 2.5;
        molecularSpeedRatioPoints[ 2 ] = 5.0;
        molecularSpeedRatioPoints[ 3 ] = 10.0;
    }
    // Set default points for high speed analysis.
    else if ( molecularSpeedRatioRegime == "High" )
    {
        molecularSpeedRatioPoints.resize( 4 );
        molecularSpeedRatioPoints[ 0 ] = 10.0;
        molecularSpeedRatioPoints[ 1 ] = 25.0;
        molecularSpeedRatioPoints[ 2 ] = 50.0;
        molecularSpeedRatioPoints[ 3 ] = 100.0;
    }
    return molecularSpeedRatioPoints;
}

//! Returns default values of angle of attack for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowAngleOfAttackPoints(
        const std::string& angleOfAttackRegime )
{
    std::vector< double > angleOfAttackPoints;

    // Set default angles of attack
    double a = - 35;
    while ( a <= 35 )
    {
        angleOfAttackPoints.push_back( convertDegreesToRadians( a ) );
        a += 5;
    }

    // Add extra points if required
    if ( angleOfAttackRegime == "Full" )
    {
        std::vector< double > frontExtension = { convertDegreesToRadians( -85.0 ),
                                                 convertDegreesToRadians( -70.0 ),
                                                 convertDegreesToRadians( -55.0 ),
                                                 convertDegreesToRadians( -40.0 ) };
        std::vector< double > rearExtension = { convertDegreesToRadians( 40.0 ),
                                                convertDegreesToRadians( 55.0 ),
                                                convertDegreesToRadians( 70.0 ),
                                                convertDegreesToRadians( 85.0 ) };
        angleOfAttackPoints.insert( angleOfAttackPoints.begin( ), frontExtension.begin( ), frontExtension.end( ) );
        angleOfAttackPoints.insert( angleOfAttackPoints.end( ), rearExtension.begin( ), rearExtension.end( ) );
    }
    return angleOfAttackPoints;
}

//! Default constructor.
RarefiedFlowAnalysis::RarefiedFlowAnalysis(
        const std::string& SPARTAExecutable,
        const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
        boost::shared_ptr< TabulatedAtmosphere > atmosphereModel,
        const std::string& simulationGases,
        const std::string& geometryFileUser,
        const double referenceArea,
        const double referenceLength,
        const int referenceAxis,
        const Eigen::Vector3d& momentReferencePoint,
        const double gridSpacing,
        const double simulatedParticlesPerCell,
        const double wallTemperature,
        const double accommodationCoefficient )
    : AerodynamicCoefficientGenerator< 3, 6 >(
          dataPointsOfIndependentVariables, referenceLength, referenceArea, referenceLength,
          momentReferencePoint,
          boost::assign::list_of( number_density_dependent )( molecular_speed_ratio_dependent )( angle_of_attack_dependent ),
          true, true ),
      SPARTAExecutable_( SPARTAExecutable ),simulationGases_( simulationGases ), referenceAxis_( referenceAxis ),
      gridSpacing_( gridSpacing ), simulatedParticlesPerCell_( simulatedParticlesPerCell ),
      wallTemperature_( wallTemperature ), accommodationCoefficient_( accommodationCoefficient )
{
    // Analyze vehicle geometry
    analyzeGeometryFile( geometryFileUser );

    // Find atmospheric conditions based on number density
    for ( unsigned int n = 0; n < dataPointsOfIndependentVariables_.at( 0 ).size( ); n++ )
    {
        atmosphericConditions_[ density_index ].push_back( atmosphereModel->getDensity( dataPointsOfIndependentVariables_.at( 0 ).at( n ) ) );
        atmosphericConditions_[ pressure_index ].push_back( atmosphereModel->getPressure( dataPointsOfIndependentVariables_.at( 0 ).at( n ) ) );
        atmosphericConditions_[ temperature_index ].push_back( atmosphereModel->getTemperature( dataPointsOfIndependentVariables_.at( 0 ).at( n ) ) );
        atmosphericConditions_[ speed_of_sound_index ].push_back( atmosphereModel->getSpeedOfSound( dataPointsOfIndependentVariables_.at( 0 ).at( n ) ) );
        atmosphericConditions_[ number_density_index ].push_back( tudat::physical_constants::AVOGADRO_CONSTANT / tudat::physical_constants::MOLAR_GAS_CONSTANT *
                                                                  atmosphericConditions_[ density_index ].at( n ) *
                                                                  atmosphereModel->getSpecificGasConstant( dataPointsOfIndependentVariables_.at( 0 ).at( n ) ) );
    }

    // Get simulation conditions
    getSimulationConditions( );

    // Read SPARTA input template
    inputTemplate_ = readSpartaInputFileTemplate( );

    // Copy input shape file to default name
    std::string commandString = "cp " + geometryFileUser + " " + getSpartaInternalGeometryFile( );
    std::system( commandString.c_str( ) );

    // Run SPARTA simulation
    generateCoefficients( );

    // Create interpolator object
    createInterpolator( );
}

//! Get aerodynamic coefficients.
void RarefiedFlowAnalysis::analyzeGeometryFile( const std::string& geometryFileUser )
{
    // Extract information on vehicle geometry
    std::pair< Eigen::Matrix< double, Eigen::Dynamic, 3 >, Eigen::Matrix< int, Eigen::Dynamic, 3 > >
            geometryData = readSpartaGeometryFile( geometryFileUser );
    shapePoints_ = geometryData.first;
    shapeTriangles_ = geometryData.second;
    numberOfPoints_ = shapePoints_.rows( );
    numberOfTriangles_ = shapeTriangles_.rows( );

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
        elementMomentArm_.col( i ) = currentCentroid - momentReferencePoint_;
    }

    // Compute cross-sectional area
    for ( unsigned int i = 0; i < 3; i++ )
    {
        shapeCrossSectionalArea_( i ) =
                0.5 * elementSurfaceNormal_.row( i ).cwiseAbs( ).dot( elementSurfaceArea_ );
    }

    // Check consistency with input dimensions
    const double tolerance = 1e-5;
    if ( std::fabs( shapeCrossSectionalArea_( static_cast< unsigned int >( referenceAxis_ ) ) -
                    referenceArea_ ) > tolerance )
    {
        std::cout << shapeCrossSectionalArea_( static_cast< unsigned int >( referenceAxis_ ) ) -
                     referenceArea_ << std::endl;
        throw std::runtime_error( "Error in SPARTA geometry file. Input reference area does not match the combination of "
                                  "reference axis and geometry. Tolerance set to: " + std::to_string( tolerance ) );
    }
}

//! Generate aerodynamic database.
void RarefiedFlowAnalysis::getSimulationConditions( )
{
    // Simulation boundary and grid
    for ( unsigned int i = 0; i < 3; i++ )
    {
        simulationBoundaries_( 2 * i ) = 1.5 * minimumDimensions_( i ); // add extra space around shape
        simulationBoundaries_( 2 * i + 1 ) = 1.5 * maximumDimensions_( i ); // add extra space around shape
        if ( i == static_cast< unsigned int >( referenceAxis_ ) )
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
    for ( unsigned int n = 0; n < dataPointsOfIndependentVariables_.at( 0 ).size( ); n++ )
    {
        for ( unsigned int s = 0; s < dataPointsOfIndependentVariables_.at( 1 ).size( ); s++ )
        {
            freeStreamVelocities_( n, s ) = dataPointsOfIndependentVariables_.at( 1 ).at( s ) *
                    atmosphericConditions_[ speed_of_sound_index ].at( n );
            simulationTimeStep_( n, s ) = 0.1 * ( maximumDimensions_( static_cast< unsigned int >( referenceAxis_ ) ) -
                                                  minimumDimensions_( static_cast< unsigned int >( referenceAxis_ ) ) ) /
                    freeStreamVelocities_( n, s );
            // time step is taken as time it takes for a particle to travel for 10 % of the box
        }
        ratioOfRealToSimulatedParticles_( n ) = atmosphericConditions_[ number_density_index ].at( n ) *
                std::pow( gridSpacing_, 3 ) / simulatedParticlesPerCell_;
    }
}

//! Generate aerodynamic coefficients at a single set of independent variables.
void RarefiedFlowAnalysis::generateCoefficients( )
{
    // Generate command string for SPARTA
    std::cout << "Initiating SPARTA simulation. This may take a while." << std::endl;
    std::string runSPARTACommandString = "cd " + getSpartaDataPath( ) + "; " +
            SPARTAExecutable_ + " -echo log -screen none -in " + getSpartaInputFile( );

    // Predefine variables
    Eigen::Vector3d velocityVector;
    std::string temporaryOutputFile = getSpartaOutputPath( ) + "/" + "coeff";
    std::vector< std::string > outputFileExtensions = { ".400", ".600", ".800", ".1000" };
    Eigen::Matrix< double, Eigen::Dynamic, 7 > outputMatrix;
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanPressureValues;
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanShearValues;

    // Loop over simulation parameters and run SPARTA
    // Loop over number density
    for ( unsigned int n = 0; n < dataPointsOfIndependentVariables_.at( 0 ).size( ); n++ )
    {
        // Show progress
        std::cout << std::endl << "Number density: "
                  << dataPointsOfIndependentVariables_.at( 0 ).at( n ) / 1e3
                  << std::endl;

        // Loop over molecular speed ratios
        for ( unsigned int s = 0; s < dataPointsOfIndependentVariables_.at( 1 ).size( ); s++ )
        {
            // Show progress
            std::cout << "Molecular speed ratio: "
                      << dataPointsOfIndependentVariables_.at( 1 ).at( s )
                      << std::endl;

            // Loop over angles of attack
            for ( unsigned int a = 0; a < dataPointsOfIndependentVariables_.at( 2 ).size( ); a++ )
            {
                // Show progress
                std::cout << "Angle of attack: "
                          << convertRadiansToDegrees( dataPointsOfIndependentVariables_.at( 2 ).at( a ) )
                          << " deg" << std::endl;

                // Get velocity vector
                velocityVector = Eigen::Vector3d::Zero( );
                velocityVector( static_cast< unsigned int >( referenceAxis_ ) ) = ( std::signbit( referenceAxis_ ) ? 1.0 : -1.0 ) *
                        freeStreamVelocities_( n, s );

                // Print to file
                FILE * fileIdentifier = std::fopen( getSpartaInputFile( ).c_str( ), "w" );
                std::fprintf( fileIdentifier, inputTemplate_.c_str( ),
                              simulationBoundaries_( 0 ), simulationBoundaries_( 1 ), simulationBoundaries_( 2 ),
                              simulationBoundaries_( 3 ), simulationBoundaries_( 4 ), simulationBoundaries_( 5 ),
                              simulationGrid_( 0 ), simulationGrid_( 1 ), simulationGrid_( 2 ),
                              atmosphericConditions_[ number_density_index ].at( n ), ratioOfRealToSimulatedParticles_( n ),
                              simulationGases_.c_str( ),
                              simulationGases_.c_str( ), velocityVector( 0 ), velocityVector( 1 ), velocityVector( 2 ),
                              simulationGases_.c_str( ), atmosphericConditions_[ temperature_index ].at( n ),
                              convertRadiansToDegrees( dataPointsOfIndependentVariables_.at( 2 ).at( a ) ),
                              wallTemperature_, accommodationCoefficient_,
                              simulationTimeStep_( n, s ),
                              getSpartaOutputDirectory( ).c_str( ) );
                std::fclose( fileIdentifier );

                // Run SPARTA
                int systemStatus = std::system( runSPARTACommandString.c_str( ) );
                if ( systemStatus != 0 )
                {
                    throw std::runtime_error( "Error: SPARTA simulation failed. See the log.sparta file in "
                                              "Tudat/External/SPARTA/ for more details." );
                }

                // Loop over angles of attack
                meanPressureValues.resize( 3, numberOfTriangles_ );
                meanShearValues.resize( 3, numberOfTriangles_ );

                // Read output files and compute mean pressure and shear force values
                meanPressureValues.setZero( );
                meanShearValues.setZero( );
                for ( unsigned int i = 0; i < outputFileExtensions.size( ); i++ )
                {
                    outputMatrix = readMatrixFromFile( temporaryOutputFile + outputFileExtensions.at( i ), "\t ;,", "%", 9 );
                    for ( unsigned int j = 0; j < 3; j++ )
                    {
                        meanPressureValues.row( j ) += outputMatrix.col( j + 1 ).transpose( );
                        meanShearValues.row( j ) += outputMatrix.col( j + 4 ).transpose( );
                    }
                }
                meanPressureValues /= outputFileExtensions.size( );
                meanShearValues /= outputFileExtensions.size( );

                // Convert pressure and shear forces to aerodynamic coefficients
                aerodynamicCoefficients_[ n ][ s ][ a ] = computeAerodynamicCoefficientsFromPressureShearForces(
                            meanPressureValues,
                            meanShearValues,
                            atmosphericConditions_[ density_index ].at( n ),
                            atmosphericConditions_[ pressure_index ].at( n ),
                            freeStreamVelocities_( n, s ),
                            elementSurfaceNormal_,
                            elementSurfaceArea_,
                            elementMomentArm_,
                            referenceArea_,
                            referenceLength_ );

                // Clean up results folder
                std::string commandString = "rm " + getSpartaOutputPath( ) + "/coeff.*";
                std::system( commandString.c_str( ) );
            }
        }
    }
}

//! Get aerodynamic coefficients at specific conditions.
Eigen::Vector6d RarefiedFlowAnalysis::getAerodynamicCoefficientsDataPoint(
        const boost::array< int, 3 > independentVariables )
{
    // Return requested coefficients.
    return aerodynamicCoefficients_( independentVariables );
}

} // namespace aerodynamics

} // namespace tudat
