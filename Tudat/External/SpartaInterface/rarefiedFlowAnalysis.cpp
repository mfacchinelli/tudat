/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Geometry>

#include "Tudat/External/SpartaInterface/rarefiedFlowAnalysis.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

namespace tudat
{

namespace aerodynamics
{

using namespace unit_conversions;

//! Returns default values of altitude for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowAltitudePoints( const std::string& targetPlanet )
{
    std::vector< double > altitudePoints;

    // Set default points for Earth.
    if ( targetPlanet == "Earth" )
    {
        altitudePoints.resize( 5 );
        altitudePoints[ 0 ] = 225.0e3;
        altitudePoints[ 1 ] = 250.0e3;
        altitudePoints[ 2 ] = 300.0e3;
        altitudePoints[ 3 ] = 400.0e3;
        altitudePoints[ 4 ] = 600.0;
    }
    // Set default points for Mars.
    else if ( targetPlanet == "Mars" )
    {
        altitudePoints.resize( 5 );
        altitudePoints[ 0 ] = 125.0e3;
        altitudePoints[ 1 ] = 150.0e3;
        altitudePoints[ 2 ] = 200.0e3;
        altitudePoints[ 3 ] = 300.0e3;
        altitudePoints[ 4 ] = 500.0e3;
    }
    // Give error otherwise.
    else
    {
        throw std::runtime_error( "Error in altitude range selection for SPARTA simulation. "
                                  "Planet not supported. Use a custom altitude range instead." );
    }
    return altitudePoints;
}

//! Returns default values of Mach number for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowMachPoints( const std::string& machRegime )
{
    std::vector< double > machPoints;

    // Set default points for full hypersonic analysis.
    if ( machRegime == "Full" )
    {
        machPoints.resize( 6 );

        machPoints[ 0 ] = 3.0;
        machPoints[ 1 ] = 4.0;
        machPoints[ 2 ] = 5.0;
        machPoints[ 3 ] = 8.0;
        machPoints[ 4 ] = 10.0;
        machPoints[ 5 ] = 20.0;
    }
    // Set default points for low hypersonic analysis.
    else if ( machRegime == "Low" )
    {
        machPoints.resize( 5 );
        machPoints[ 0 ] = 3.0;
        machPoints[ 1 ] = 4.0;
        machPoints[ 2 ] = 5.0;
        machPoints[ 3 ] = 8.0;
        machPoints[ 4 ] = 10.0;
    }
    // Set default points for high hypersonic analysis.
    else if ( machRegime == "High" )
    {
        machPoints.resize( 4 );
        machPoints[ 0 ] = 5.0;
        machPoints[ 1 ] = 8.0;
        machPoints[ 2 ] = 10.0;
        machPoints[ 3 ] = 20.0;
    }
    return machPoints;
}

//! Returns default values of angle of attack for use in RarefiedFlowAnalysis.
std::vector< double > getDefaultRarefiedFlowAngleOfAttackPoints( const std::string& angleOfAttackRegime )
{
    std::vector< double > angleOfAttackPoints;

    // Set default angles of attack
    double a = -35.0;
    while ( a <= 35.0 )
    {
        angleOfAttackPoints.push_back( convertDegreesToRadians( a ) );
        a += 5.0;
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
        const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
        boost::shared_ptr< TabulatedAtmosphere > atmosphereModel,
        const std::string& simulationGases, const std::string& geometryFileUser,
        const double referenceArea, const double referenceLength, const int referenceAxis,
        const Eigen::Vector3d& momentReferencePoint,
        const double gridSpacing, const double simulatedParticlesPerCell,
        const double wallTemperature, const double accommodationCoefficient,
        const bool printProgressInCommandWindow,
        const std::string& SpartaExecutable, const std::string& MpiExecutable, const unsigned int numberOfCores ) :
    AerodynamicCoefficientGenerator< 3, 6 >(
        dataPointsOfIndependentVariables, referenceLength, referenceArea, referenceLength,
        momentReferencePoint,
        boost::assign::list_of( altitude_dependent )( mach_number_dependent )( angle_of_attack_dependent ),
        true, true ),
    sparta_interface::SpartaInterface( dataPointsOfIndependentVariables, simulationGases, referenceAxis,
                                       gridSpacing, simulatedParticlesPerCell,
                                       wallTemperature, accommodationCoefficient, printProgressInCommandWindow,
                                       SpartaExecutable, MpiExecutable, numberOfCores )
{
    // Check that reference dimension makes sense
    if ( referenceDimension_ > 2 )
    {
        throw std::runtime_error( "Error in SPARTA rarefied flow analysis. Reference axis makes "
                                  "no sense for a universe with 3 spacial dimensions (i.e., our universe). "
                                  "Note that the first dimension is identified with 0." );
    }

    // Analyze vehicle geometry
    this->analyzeGeometryFile( geometryFileUser, referenceArea_, momentReferencePoint_ );

    // Retrieve atmospheric conditions based on altitude
    for ( unsigned int h = 0; h < dataPointsOfIndependentVariables_.at( 0 ).size( ); h++ )
    {
        atmosphericConditions_[ sparta_interface::density_index ].push_back(
                    atmosphereModel->getDensity( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ sparta_interface::pressure_index ].push_back(
                    atmosphereModel->getPressure( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ sparta_interface::temperature_index ].push_back(
                    atmosphereModel->getTemperature( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ sparta_interface::speed_of_sound_index ].push_back(
                    atmosphereModel->getSpeedOfSound( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
        atmosphericConditions_[ sparta_interface::number_density_index ].push_back(
                    tudat::physical_constants::AVOGADRO_CONSTANT * atmosphericConditions_[ sparta_interface::density_index ].at( h ) /
                    atmosphereModel->getMolarMass( dataPointsOfIndependentVariables_.at( 0 ).at( h ) ) );
    }

    // Get simulation conditions
    this->getSimulationConditions( );

    // Run SPARTA simulation
    generateCoefficients( );

    // Create interpolator object
    this->createInterpolator( );
}

//! Generate aerodynamic database.
void RarefiedFlowAnalysis::generateCoefficients( )
{
    // Inform user on progress
    std::cout << "Initiating SPARTA simulation. This may take a while." << std::endl;

    // Loop over simulation parameters and run SPARTA
    // Loop over altitude
    for ( unsigned int h = 0; h < dataPointsOfIndependentVariables_.at( 0 ).size( ); h++ )
    {
        // Inform user on progress
        std::cout << std::endl << "Altitude: "
                  << dataPointsOfIndependentVariables_.at( 0 ).at( h ) / 1e3
                  << " km" << std::endl;

        // Loop over Mach numbers
        for ( unsigned int m = 0; m < dataPointsOfIndependentVariables_.at( 1 ).size( ); m++ )
        {
            // Inform user on progress
            std::cout << "Mach number: "
                      << dataPointsOfIndependentVariables_.at( 1 ).at( m )
                      << std::endl;

            // Loop over angles of attack
            for ( unsigned int a = 0; a < dataPointsOfIndependentVariables_.at( 2 ).size( ); a++ )
            {
                // Inform user on progress
                std::cout << "Angle of attack: "
                          << convertRadiansToDegrees( dataPointsOfIndependentVariables_.at( 2 ).at( a ) )
                          << " deg" << std::endl;

                // Run SPARTA simulation
                this->runSpartaSimulation( h, m, a );

                // Process SPARTA output
                aerodynamicCoefficients_[ h ][ m ][ a ] = this->processSpartaOutput( referenceArea_, referenceLength_, h, m );
            }
        }
    }

    // Inform user on progress
    std::cout << std::endl << "SPARTA simulation complete." << std::endl << std::endl;
}

//! Get aerodynamic coefficients at specific conditions.
Eigen::Vector6d RarefiedFlowAnalysis::getAerodynamicCoefficientsDataPoint( const boost::array< int, 3 > independentVariables )
{
    // Return requested coefficients.
    return aerodynamicCoefficients_( independentVariables );
}

} // namespace aerodynamics

} // namespace tudat
