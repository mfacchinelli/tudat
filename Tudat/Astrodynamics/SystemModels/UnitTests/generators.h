/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Facchinelli, M. (2018). Aerobraking Navigation, Guidance and Control.
 *          Master Thesis, Delft University of Technology.
 */

#ifndef TUDAT_GENERATORS_H
#define TUDAT_GENERATORS_H

#include <limits>
#include <vector>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

#include "Tudat/Astrodynamics/SystemModels/controlSystem.h"
#include "Tudat/Astrodynamics/SystemModels/guidanceSystem.h"
#include "Tudat/Astrodynamics/SystemModels/instrumentsModel.h"
#include "Tudat/Astrodynamics/SystemModels/navigationSystem.h"

namespace tudat
{

namespace system_models
{

//! Class to generate instruments model.
class InstrumentsModelGenerator
{
public:

    //! Constructor
    InstrumentsModelGenerator( )
    {
        // Load Spice kernels
        spice_interface::loadStandardSpiceKernels( );

        // Simulation times
        const double simulationStartEpoch = 7.0 * tudat::physical_constants::JULIAN_YEAR +
                30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
        const double simulationEndEpoch = simulationStartEpoch + 1.0;

        // Set body names
        spacecraftName_ = "Satellite";
        planetName_ = "Mars";
        std::string sunBody = "Sun";
        std::string thirdBody = "Earth";
        std::vector< std::string > bodiesToCreate = { planetName_, sunBody, thirdBody };

        ///////////////////////////////////////////////  Create body map   ///////////////////////////////////////////////////////////////

        // Create body settings and give default values
        std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettings =
                simulation_setup::getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
        for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        }
        bodySettings[ thirdBody ]->gravityFieldSettings = boost::make_shared< simulation_setup::CentralGravityFieldSettings >( 0.0 );

        // Create body objects
        bodyMap_ = createBodies( bodySettings );

        // Create spacecraft object
        bodyMap_[ spacecraftName_ ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap_[ spacecraftName_ ]->setConstantBodyMass( 1000.0 );

        // Set radiation pressure settings
        boost::shared_ptr< simulation_setup::RadiationPressureInterfaceSettings > radiationPressureSettings =
                boost::make_shared< simulation_setup::CannonBallRadiationPressureInterfaceSettings >( sunBody, 10.0, 1.0 );
        bodyMap_[ spacecraftName_ ]->setRadiationPressureInterface( sunBody, createRadiationPressureInterface(
                                                                        radiationPressureSettings, spacecraftName_, bodyMap_ ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap_, "SSB", "J2000" );

        ///////////////////////////////////////////////  Create acceleration map   ///////////////////////////////////////////////////////

        // Define acceleration settings for simulation model
        std::map< std::string, std::vector< boost::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::SphericalHarmonicAccelerationSettings >( 21, 21 ) );
        for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            if ( bodiesToCreate.at( i ) != planetName_ )
            {
                accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back(
                            boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::central_gravity ) );
            }
        }
        accelerationsOfSatellite[ sunBody ].push_back(
                    boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );

        // Set accelerations settings
        simulation_setup::SelectedAccelerationMap accelerationMap;
        accelerationMap[ spacecraftName_ ] = accelerationsOfSatellite;
        accelerationModelMap_ = simulation_setup::createAccelerationModelsMap( bodyMap_, accelerationMap, { spacecraftName_ }, { planetName_ } );

        ///////////////////////////////////////////////  Create instrument model   ///////////////////////////////////////////////////

        // Declare instrument model
        instrumentsModel_ = boost::make_shared< InstrumentsModel >( bodyMap_, accelerationModelMap_, spacecraftName_, planetName_ );

        ///////////////////////////////////////////////  Create simulation dynamics   ////////////////////////////////////////////////

        // Set initial conditions
        Eigen::Vector6d initialSpacecraftKeplerianState;
        initialSpacecraftKeplerianState[ 0 ] = 10000.0e3;
        initialSpacecraftKeplerianState[ 1 ] = 0.1;
        initialSpacecraftKeplerianState[ 2 ] = unit_conversions::convertDegreesToRadians( 37.6 );
        initialSpacecraftKeplerianState[ 3 ] = 0.0;
        initialSpacecraftKeplerianState[ 4 ] = 0.0;
        initialSpacecraftKeplerianState[ 5 ] = 0.0;
        const Eigen::Vector6d initialSpacecraftCartesianState = orbital_element_conversions::convertKeplerianToCartesianElements(
                    initialSpacecraftKeplerianState, bodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( ) );

        // Create integrator settings
        integratorSettings_ = boost::make_shared< numerical_integrators::IntegratorSettings< > >(
                    numerical_integrators::rungeKutta4, simulationStartEpoch, 1.0 );

        // Create propagator settings
        propagatorSettings_ = boost::make_shared< propagators::TranslationalStatePropagatorSettings< > >(
                    std::vector< std::string >{ planetName_ }, accelerationModelMap_,
                    std::vector< std::string >{ spacecraftName_ }, initialSpacecraftCartesianState,
                    boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ),
                    propagators::cowell );

        // Create empty dynamics simulator object
        dynamicsSimulator_ = boost::make_shared< propagators::SingleArcDynamicsSimulator< > >(
                    bodyMap_, integratorSettings_, propagatorSettings_, false );
    }

    //! Destructor
    ~InstrumentsModelGenerator( ) { }

    //! Function to update the environment.
    void updateEnvironment( )
    {
        // Propagate state and retrieve solution
        dynamicsSimulator_->integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );
        std::map< double, Eigen::VectorXd > numericalSolution = dynamicsSimulator_->getEquationsOfMotionNumericalSolution( );

        // Reset integrator and propagator settings
        integratorSettings_->initialTime_ = numericalSolution.rbegin( )->first;
        propagatorSettings_->resetInitialStates( numericalSolution.rbegin( )->second );
        boost::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< > >(
                    propagatorSettings_ )->resetTerminationSettings(
                    boost::make_shared< propagators::PropagationTimeTerminationSettings >( integratorSettings_->initialTime_ +
                                                                                           integratorSettings_->initialTimeStep_ ) );

        // Reset dynamics simulator object
        dynamicsSimulator_ = boost::make_shared< propagators::SingleArcDynamicsSimulator< > >(
                    bodyMap_, integratorSettings_, propagatorSettings_, false );

        // Update onboard models
        instrumentsModel_->updateInstruments( integratorSettings_->initialTime_ );
    }

    //! Function to retrive instruments model.
    boost::shared_ptr< InstrumentsModel > getInstrumentsModel( )
    {
        return instrumentsModel_;
    }

private:

    //! Name of spacecraft.
    std::string spacecraftName_;

    //! Name of planet.
    std::string planetName_;

    //! Body map pointer.
    tudat::simulation_setup::NamedBodyMap bodyMap_;

    //! Acceleration model map pointer.
    tudat::basic_astrodynamics::AccelerationMap accelerationModelMap_;

    //! Instruments model pointer.
    boost::shared_ptr< InstrumentsModel > instrumentsModel_;

    //! Dynamics simulator pointer.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > >  integratorSettings_;

    //! Dynamics simulator pointer.
    boost::shared_ptr< propagators::PropagatorSettings< > > propagatorSettings_;

    //! Dynamics simulator pointer.
    boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator_;

};

//! Class to generate navigation system.
class NavigationSystemGenerator
{
public:

    //! Constructor
    NavigationSystemGenerator( const aerodynamics::AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel ) :
        selectedOnboardAtmosphereModel_( selectedOnboardAtmosphereModel )
    {
        // Load Spice kernels
        spice_interface::loadStandardSpiceKernels( );

        // Simulation times
        simulationStartEpoch_ = 7.0 * tudat::physical_constants::JULIAN_YEAR +
                30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
        const double simulationEndEpoch = simulationStartEpoch_ + 1.4 * tudat::physical_constants::JULIAN_DAY;

        // Set body names
        spacecraftName_ = "Satellite";
        planetName_ = "Mars";

        ///////////////////////////////////////////////  Create onboard body map   ///////////////////////////////////////////////////

        // Create body settings and give default values
        std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > onboardBodySettings =
                simulation_setup::getDefaultBodySettings( { planetName_ }, simulationStartEpoch_ - 300.0, simulationEndEpoch + 300.0 );
        onboardBodySettings[ planetName_ ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        onboardBodySettings[ planetName_ ]->rotationModelSettings->resetOriginalFrame( "J2000" );

        std::vector< double > vectorOfModelSpecificParameters;
        switch ( selectedOnboardAtmosphereModel_ )
        {
        case aerodynamics::exponential_atmosphere_model:
            vectorOfModelSpecificParameters = { 115.0e3, 2.424e-08, 6533.0 };
            break;
        case aerodynamics::three_term_atmosphere_model:
            vectorOfModelSpecificParameters = { 115.0e3, 2.424e-08, 6533.0, -1.0, 0.0, 0.0 };
            break;
        default:
            throw std::runtime_error( "Error in simulation. Selected atmosphere model not supported." );
        }
        onboardBodySettings[ planetName_ ]->atmosphereSettings =
                boost::make_shared< simulation_setup::CustomConstantTemperatureAtmosphereSettings >(
                    selectedOnboardAtmosphereModel_, 215.0, 197.0, 1.3, vectorOfModelSpecificParameters );

        // Create body objects
        onboardBodyMap_ = createBodies( onboardBodySettings );

        // Spacecraft parameters
        const double referenceAreaAerodynamic = 37.5;
        Eigen::Vector3d onboardAerodynamicCoefficients = Eigen::Vector3d::Zero( );
        onboardAerodynamicCoefficients[ 0 ] = 1.9;

        // Create spacecraft object
        onboardBodyMap_[ spacecraftName_ ] = boost::make_shared< simulation_setup::Body >( );
        onboardBodyMap_[ spacecraftName_ ]->setConstantBodyMass( 1000.0 );
        onboardBodyMap_[ spacecraftName_ ]->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( boost::make_shared< simulation_setup::ConstantAerodynamicCoefficientSettings >(
                                                               referenceAreaAerodynamic, onboardAerodynamicCoefficients, true, true ),
                                                           spacecraftName_ ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( onboardBodyMap_, planetName_, "J2000" );

        ///////////////////////////////////////////////  Create onboard acceleration map   ///////////////////////////////////////////

        // Define acceleration settings for simulation model
        std::map< std::string, std::vector< boost::shared_ptr< simulation_setup::AccelerationSettings > > > onboardAccelerationsOfSatellite;
        onboardAccelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::SphericalHarmonicAccelerationSettings >( 4, 4 ) );
        onboardAccelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::aerodynamic ) );

        // Set accelerations settings
        simulation_setup::SelectedAccelerationMap onboardAccelerationMap;
        onboardAccelerationMap[ spacecraftName_ ] = onboardAccelerationsOfSatellite;
        onboardAccelerationModelMap_ = simulation_setup::createAccelerationModelsMap(
                    onboardBodyMap_, onboardAccelerationMap, { spacecraftName_ }, { planetName_ } );
    }

    //! Function to generate a navigation system.
    boost::shared_ptr< NavigationSystem > createNavigationSystem(
            const Eigen::Vector9d initialEstimatedStateVector,
            const double onboardComputerRefreshStepSize )
    {
        ///////////////////////////////////////////////  Create navigation filter   //////////////////////////////////////////////////////

        // Set navigation variables
        Eigen::Matrix9d initialEstimatedStateCovarianceMatrix = Eigen::Matrix9d::Identity( );

        const double positionStandardDeviation = 1.0e2;
        const double translationalVelocityStandardDeviation = 1.0e-1;
        const double accelerometerBiasStandardDeviation = 1.0e-4;
        const Eigen::Vector3d positionAccuracy = Eigen::Vector3d::Constant( 5.0e2 );

        Eigen::Vector9d diagonalOfSystemUncertainty;
        diagonalOfSystemUncertainty << Eigen::Vector3d::Constant( std::pow( positionStandardDeviation, 2 ) ),
                Eigen::Vector3d::Constant( std::pow( translationalVelocityStandardDeviation, 2 ) ),
                Eigen::Vector3d::Constant( std::pow( accelerometerBiasStandardDeviation, 2 ) );
        Eigen::Matrix9d systemUncertainty = diagonalOfSystemUncertainty.asDiagonal( );

        Eigen::Matrix3d measurementUncertainty = 10.0 * ( positionAccuracy.cwiseProduct( positionAccuracy ) ).asDiagonal( );

        // Create navigation filter integration settings
        boost::shared_ptr< numerical_integrators::IntegratorSettings< > > filterIntegratorSettings =
                boost::make_shared< numerical_integrators::IntegratorSettings< > >(
                    numerical_integrators::euler, simulationStartEpoch_, onboardComputerRefreshStepSize );

        // Create navigation filter
        boost::shared_ptr< filters::FilterSettings< > > navigationFilterSettings = boost::make_shared<
                filters::UnscentedKalmanFilterSettings< > >(
                    systemUncertainty, measurementUncertainty, onboardComputerRefreshStepSize, simulationStartEpoch_,
                    initialEstimatedStateVector, initialEstimatedStateCovarianceMatrix, filterIntegratorSettings );

        ///////////////////////////////////////////////  Create navigation system   //////////////////////////////////////////////////////

        // Other variables
        const double atmosphericInterfaceAltitude = 200.0e3;
        const double reducedAtmosphericInterfaceAltitude = 150.0e3;
        const double periapseEstimatorConstant = 0.955;
        const unsigned int frequencyOfDeepSpaceNetworkTracking = -1;
        const unsigned int numberOfRequiredAtmosphereSamplesForInitiation = 7;

        // Create navigation system
        navigationSystem_ = boost::make_shared< NavigationSystem >(
                    onboardBodyMap_, onboardAccelerationModelMap_, spacecraftName_, planetName_,
                    navigationFilterSettings, selectedOnboardAtmosphereModel_,
                    atmosphericInterfaceAltitude, reducedAtmosphericInterfaceAltitude,
                    periapseEstimatorConstant, numberOfRequiredAtmosphereSamplesForInitiation,
                    frequencyOfDeepSpaceNetworkTracking, std::vector< Eigen::Vector3d >( ), reference_frames::inertial_frame,
                    std::make_pair( TUDAT_NAN, TUDAT_NAN ), true );

        // Complete navigation system
        navigationSystem_->createNavigationSystemObjects( 1, boost::lambda::constant( Eigen::Vector3d::Zero( ) ) );

        // Give ouput
        return navigationSystem_;
    }

private:

    //! Onboard atmosphere model.
    aerodynamics::AvailableConstantTemperatureAtmosphereModels selectedOnboardAtmosphereModel_;

    //! Simulation start epoch.
    double simulationStartEpoch_;

    //! Name of spacecraft.
    std::string spacecraftName_;

    //! Name of planet.
    std::string planetName_;

    //! Body map pointer.
    tudat::simulation_setup::NamedBodyMap onboardBodyMap_;

    //! Acceleration model map pointer.
    tudat::basic_astrodynamics::AccelerationMap onboardAccelerationModelMap_;

    //! Navigation system pointer.
    boost::shared_ptr< NavigationSystem > navigationSystem_;

};

//! Function to propagate state for one orbit in an aerobraking scenario.
std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagateStateForAerobrakingScenario(
        const unsigned int initialConditions, const bool useSunGravity = true )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::ephemerides;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::gravitation;
    using namespace tudat::input_output;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::simulation_setup;
    using namespace tudat::system_models;
    using namespace tudat::unit_conversions;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings
    const double simulationStartEpoch = 7.0 * tudat::physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
    double simulationEndEpoch = simulationStartEpoch;
    switch ( initialConditions )
    {
    case 0:
        simulationEndEpoch += 1.4 * tudat::physical_constants::JULIAN_DAY;
        break;
    case 1:
        simulationEndEpoch += 0.18 * tudat::physical_constants::JULIAN_DAY;
        break;
    case 2:
        simulationEndEpoch += 0.113 * tudat::physical_constants::JULIAN_DAY;
        break;
    }

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Earth" );

    // Tabulated atmosphere settings
    std::map< int, std::string > tabulatedAtmosphereFiles;
    tabulatedAtmosphereFiles[ 0 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";
    tabulatedAtmosphereFiles[ 1 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
    tabulatedAtmosphereFiles[ 2 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
    tabulatedAtmosphereFiles[ 3 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
    tabulatedAtmosphereFiles[ 4 ] = getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
    std::vector< AtmosphereIndependentVariables > atmosphereIndependentVariables = {
        longitude_dependent_atmosphere, latitude_dependent_atmosphere, altitude_dependent_atmosphere };
    std::vector< AtmosphereDependentVariables > atmosphereDependentVariables = {
        density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere,
        gas_constant_dependent_atmosphere, specific_heat_ratio_dependent_atmosphere };
    std::vector< interpolators::BoundaryInterpolationType > boundaryConditions =
            std::vector< interpolators::BoundaryInterpolationType >( 3, interpolators::use_boundary_value );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    // Give Mars a more detailed environment
    bodySettings[ "Mars" ]->gravityFieldSettings = boost::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d );
    bodySettings[ "Mars" ]->atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables, boundaryConditions );

    // Give Earth zero gravity field such that ephemeris is created, but no acceleration
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( 0.0 );
    if ( !useSunGravity )
    {
        bodySettings[ "Sun" ]->gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >( 0.0 );
    }

    // Create body objects
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Define simulation objects
    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Satellite" );

    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Mars" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Satellite" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Satellite" ]->setConstantBodyMass( 1000.0 );

    // Spacecraft parameters
    const double referenceAreaAerodynamic = 37.5;
    const double referenceAreaRadiation = 37.5;
    const double radiationPressureCoefficient = 1.0;

    // Aerodynamic coefficients from file
    std::map< int, std::string > aerodynamicForceCoefficientFiles;
    aerodynamicForceCoefficientFiles[ 0 ] = getTudatRootPath( ) + "External/MRODragCoefficients.txt";
    aerodynamicForceCoefficientFiles[ 1 ] = getTudatRootPath( ) + "External/MROSideCoefficients.txt";
    aerodynamicForceCoefficientFiles[ 2 ] = getTudatRootPath( ) + "External/MROLiftCoefficients.txt";

    // Create aerodynamic coefficient settings
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                aerodynamicForceCoefficientFiles, referenceAreaAerodynamic,
                std::vector< AerodynamicCoefficientsIndependentVariables >{ angle_of_attack_dependent, angle_of_sideslip_dependent,
                                                                            altitude_dependent } );

    // Constant radiation pressure variables
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Set aerodynamic coefficient and radiation pressure settings
    bodyMap[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );
    bodyMap[ "Satellite" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface(
                                                               radiationPressureSettings, "Satellite", bodyMap ) );

    // Finalize body creation
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define acceleration settings
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 21, 21 ) );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        if ( bodiesToCreate.at( i ) != "Mars" )
        {
            accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        }
    }
    accelerationsOfSatellite[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfSatellite[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );

    // Set accelerations settings
    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE INITIAL CONDITIONS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Satellite.
    Eigen::Vector6d initialSpacecraftKeplerianState;
    switch ( initialConditions )
    {
    case 0:
    {
        initialSpacecraftKeplerianState( semiMajorAxisIndex ) = 25946932.3;
        initialSpacecraftKeplerianState( eccentricityIndex ) = 0.8651912;
        break;
    }
    case 1:
    {
        initialSpacecraftKeplerianState( semiMajorAxisIndex ) = 5.29046e+06;
        initialSpacecraftKeplerianState( eccentricityIndex ) = 0.337791;
        break;
    }
    case 2:
    {
        initialSpacecraftKeplerianState( semiMajorAxisIndex ) = 4.70078e+06;
        initialSpacecraftKeplerianState( eccentricityIndex ) = 0.254259;
        break;
    }
    default:
        throw std::runtime_error( "Error. Only cases 0 and 1 are supported." );
    }
    initialSpacecraftKeplerianState( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 93.0 );
    initialSpacecraftKeplerianState( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 158.7 );
    initialSpacecraftKeplerianState( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 43.6 );
    initialSpacecraftKeplerianState( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 180.0 );

    // Convert to Cartesian coordinates
    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d translationalInitialState = convertKeplerianToCartesianElements(
                initialSpacecraftKeplerianState, marsGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Dependent variables
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                      basic_astrodynamics::aerodynamic, "Satellite", "Mars", false ) );
    dependentVariables.push_back( boost::make_shared< SingleDependentVariableSaveSettings >(
                                      local_density_dependent_variable, "Satellite", "Mars" ) );

    // Propagator settings
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, translationalInitialState, simulationEndEpoch,
                unified_state_model_exponential_map, boost::make_shared< DependentVariableSaveSettings >( dependentVariables, false ) );

    // Integrator settings
    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< > >(
                simulationStartEpoch, 10.0, RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0e-5, 1.0e5, 1e-15, 1e-15 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > numericalIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableSolution = dynamicsSimulator.getDependentVariableHistory( );

    // Give output
    return std::make_pair( numericalIntegrationResult, dependentVariableSolution );
}

//! Class to generate navigation system.
class GuidanceSystemGenerator
{
public:

    //! Constructor
    GuidanceSystemGenerator(  )
    {
        // Load Spice kernels
        spice_interface::loadStandardSpiceKernels( );

        // Simulation times
        simulationStartEpoch_ = 7.0 * tudat::physical_constants::JULIAN_YEAR +
                30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
        const double simulationEndEpoch = simulationStartEpoch_ + 1.0;
        currentTime_ = simulationStartEpoch_;

        // Set body names
        spacecraftName_ = "Satellite";
        planetName_ = "Mars";
        std::string sunBody = "Sun";
        std::string thirdBody = "Earth";
        std::vector< std::string > bodiesToCreate = { planetName_, sunBody, thirdBody };

        ///////////////////////////////////////////////  Create body map   ///////////////////////////////////////////////////////////////

        // Tabulated atmosphere settings
        std::map< int, std::string > tabulatedAtmosphereFiles;
        tabulatedAtmosphereFiles[ 0 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/density.dat";
        tabulatedAtmosphereFiles[ 1 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/pressure.dat";
        tabulatedAtmosphereFiles[ 2 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/temperature.dat";
        tabulatedAtmosphereFiles[ 3 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/gasConstant.dat";
        tabulatedAtmosphereFiles[ 4 ] = input_output::getAtmosphereTablesPath( ) + "MCDMeanAtmosphereTimeAverage/specificHeatRatio.dat";
        std::vector< aerodynamics::AtmosphereIndependentVariables > atmosphereIndependentVariables = {
            aerodynamics::longitude_dependent_atmosphere, aerodynamics::latitude_dependent_atmosphere,
            aerodynamics::altitude_dependent_atmosphere };
        std::vector< aerodynamics::AtmosphereDependentVariables > atmosphereDependentVariables = {
            aerodynamics::density_dependent_atmosphere, aerodynamics::pressure_dependent_atmosphere,
            aerodynamics::temperature_dependent_atmosphere, aerodynamics::gas_constant_dependent_atmosphere,
            aerodynamics::specific_heat_ratio_dependent_atmosphere };
        std::vector< interpolators::BoundaryInterpolationType > boundaryConditions =
                std::vector< interpolators::BoundaryInterpolationType >( 3, interpolators::use_boundary_value );

        // Create body settings and give default values
        std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > bodySettings =
                simulation_setup::getDefaultBodySettings( bodiesToCreate, simulationStartEpoch_ - 300.0, simulationEndEpoch + 300.0 );
        for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
        }
        bodySettings[ thirdBody ]->gravityFieldSettings = boost::make_shared< simulation_setup::CentralGravityFieldSettings >( 0.0 );

        // Give Mars a more detailed environment
        bodySettings[ planetName_ ]->gravityFieldSettings =
                boost::make_shared< simulation_setup::FromFileSphericalHarmonicsGravityFieldSettings >( simulation_setup::jgmro120d );
        bodySettings[ planetName_ ]->atmosphereSettings = boost::make_shared< simulation_setup::TabulatedAtmosphereSettings >(
                    tabulatedAtmosphereFiles, atmosphereIndependentVariables, atmosphereDependentVariables, boundaryConditions );

        // Create body objects
        bodyMap_ = createBodies( bodySettings );

        // Create spacecraft object
        const double spacecraftMass = 1000.0;
        bodyMap_[ "Satellite" ] = boost::make_shared< simulation_setup::Body >( );
        bodyMap_[ "Satellite" ]->setConstantBodyMass( spacecraftMass );

        // Spacecraft parameters
        const double referenceAreaAerodynamic = 37.5;
        const double referenceAreaRadiation = 37.5;
        const double radiationPressureCoefficient = 1.0;

        // Aerodynamic coefficients from file
        std::map< int, std::string > aerodynamicForceCoefficientFiles;
        aerodynamicForceCoefficientFiles[ 0 ] = input_output::getTudatRootPath( ) + "External/MRODragCoefficients.txt";
        aerodynamicForceCoefficientFiles[ 1 ] = input_output::getTudatRootPath( ) + "External/MROSideCoefficients.txt";
        aerodynamicForceCoefficientFiles[ 2 ] = input_output::getTudatRootPath( ) + "External/MROLiftCoefficients.txt";

        // Create aerodynamic coefficient settings
        boost::shared_ptr< simulation_setup::AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                simulation_setup::readTabulatedAerodynamicCoefficientsFromFiles(
                    aerodynamicForceCoefficientFiles, referenceAreaAerodynamic,
                    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >{
                        aerodynamics::angle_of_attack_dependent, aerodynamics::angle_of_sideslip_dependent, aerodynamics::altitude_dependent } );

        // Constant radiation pressure variables
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( "Mars" );
        boost::shared_ptr< simulation_setup::RadiationPressureInterfaceSettings > radiationPressureSettings =
                boost::make_shared< simulation_setup::CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Set aerodynamic coefficient and radiation pressure settings
        bodyMap_[ "Satellite" ]->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );
        bodyMap_[ "Satellite" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface(
                                                                    radiationPressureSettings, "Satellite", bodyMap_ ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap_, "SSB", "J2000" );

        ///////////////////////////////////////////////  Create acceleration map   ///////////////////////////////////////////////////////

        // Define acceleration settings for simulation model
        std::map< std::string, std::vector< boost::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::SphericalHarmonicAccelerationSettings >( 21, 21 ) );
        for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            if ( bodiesToCreate.at( i ) != planetName_ )
            {
                accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back(
                            boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::central_gravity ) );
            }
        }
        accelerationsOfSatellite[ sunBody ].push_back(
                    boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );
        accelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::aerodynamic ) );

        // Set accelerations settings
        simulation_setup::SelectedAccelerationMap accelerationMap;
        accelerationMap[ spacecraftName_ ] = accelerationsOfSatellite;
        accelerationModelMap_ = simulation_setup::createAccelerationModelsMap( bodyMap_, accelerationMap, { spacecraftName_ }, { planetName_ } );

        ///////////////////////////////////////////////  Create simulation dynamics   ////////////////////////////////////////////////

        // Set initial conditions
        Eigen::Vector6d initialSpacecraftKeplerianState;
        initialSpacecraftKeplerianState[ 0 ] = 26021000.0;
        initialSpacecraftKeplerianState[ 1 ] = 0.859882;
        initialSpacecraftKeplerianState( 2 ) = unit_conversions::convertDegreesToRadians( 93.0 );
        initialSpacecraftKeplerianState( 4 ) = unit_conversions::convertDegreesToRadians( 158.7 );
        initialSpacecraftKeplerianState( 3 ) = unit_conversions::convertDegreesToRadians( 43.6 );
        initialSpacecraftKeplerianState( 5 ) = unit_conversions::convertDegreesToRadians( 180.0 );
        const Eigen::Vector6d initialSpacecraftCartesianState = orbital_element_conversions::convertKeplerianToCartesianElements(
                    initialSpacecraftKeplerianState, bodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( ) );

        // Create object of dependent variables to save
        std::vector< boost::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
        dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                              propagators::local_dynamic_pressure_dependent_variable, spacecraftName_, planetName_ ) );
        dependentVariablesList.push_back( boost::make_shared< propagators::SingleDependentVariableSaveSettings >(
                                              propagators::local_aerodynamic_heat_rate_dependent_variable, spacecraftName_, planetName_ ) );

        // Create integrator settings
        integratorSettings_ = boost::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< > >(
                    simulationStartEpoch_, 10.0, numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg78,
                    1.0e-5, 1.0e5, 1.0e-15, 1.0e-15 );

        // Create propagator settings
        propagatorSettings_ = boost::make_shared< propagators::TranslationalStatePropagatorSettings< > >(
                    std::vector< std::string >{ planetName_ }, accelerationModelMap_,
                    std::vector< std::string >{ spacecraftName_ }, initialSpacecraftCartesianState,
                    boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ),
                    propagators::unified_state_model_exponential_map,
                    boost::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList, false ) );

        ///////////////////////////////////////////////  Create onboard body map   ///////////////////////////////////////////////////

        // Create body settings and give default values
        std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > > onboardBodySettings =
                simulation_setup::getDefaultBodySettings( { planetName_ }, simulationStartEpoch_ - 300.0, simulationEndEpoch + 300.0 );
        onboardBodySettings[ planetName_ ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        onboardBodySettings[ planetName_ ]->rotationModelSettings->resetOriginalFrame( "J2000" );

        std::vector< double > vectorOfModelSpecificParameters = { 115.0e3, 2.424e-08, 6533.0 };
        onboardBodySettings[ planetName_ ]->atmosphereSettings =
                boost::make_shared< simulation_setup::CustomConstantTemperatureAtmosphereSettings >(
                    aerodynamics::exponential_atmosphere_model, 215.0, 197.0, 1.3, vectorOfModelSpecificParameters );

        // Create body objects
        onboardBodyMap_ = createBodies( onboardBodySettings );

        // Spacecraft parameters
        Eigen::Vector3d onboardAerodynamicCoefficients = Eigen::Vector3d::Zero( );
        onboardAerodynamicCoefficients[ 0 ] = 1.9;

        // Create spacecraft object
        onboardBodyMap_[ spacecraftName_ ] = boost::make_shared< simulation_setup::Body >( );
        onboardBodyMap_[ spacecraftName_ ]->setConstantBodyMass( spacecraftMass );
        onboardBodyMap_[ spacecraftName_ ]->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( boost::make_shared< simulation_setup::ConstantAerodynamicCoefficientSettings >(
                                                               referenceAreaAerodynamic, onboardAerodynamicCoefficients, true, true ),
                                                           spacecraftName_ ) );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( onboardBodyMap_, planetName_, "J2000" );

        ///////////////////////////////////////////////  Create onboard acceleration map   ///////////////////////////////////////////

        // Define acceleration settings for simulation model
        std::map< std::string, std::vector< boost::shared_ptr< simulation_setup::AccelerationSettings > > > onboardAccelerationsOfSatellite;
        onboardAccelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::SphericalHarmonicAccelerationSettings >( 4, 4 ) );
        onboardAccelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::aerodynamic ) );

        // Set accelerations settings
        simulation_setup::SelectedAccelerationMap onboardAccelerationMap;
        onboardAccelerationMap[ spacecraftName_ ] = onboardAccelerationsOfSatellite;
        onboardAccelerationModelMap_ = simulation_setup::createAccelerationModelsMap(
                    onboardBodyMap_, onboardAccelerationMap, { spacecraftName_ }, { planetName_ } );

        ///////////////////////////////////////////////  Create onboard dynamics   ///////////////////////////////////////////////////

        // Create object for propagation of spacecraft state with user-provided initial conditions
        onboardIntegratorSettings_ = boost::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    0.0, 10.0, numerical_integrators::RungeKuttaCoefficients::rungeKuttaFehlberg56, 0.1, 100.0, 1.0e-12, 1.0e-12 );
        onboardPropagatorSettings_ = boost::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ planetName_ }, onboardAccelerationModelMap_,
                    std::vector< std::string >{ spacecraftName_ }, Eigen::Vector6d::Zero( ),
                    0.0, propagators::unified_state_model_exponential_map,
                    boost::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList, false ) );
    }

    //! Function to generate a navigation system.
    boost::shared_ptr< GuidanceSystem > createGuidanceSystem( const double targetPeriapsisAltitude,
                                                              const double targetApoapsisAltitude,
                                                              const double maximumAllowedHeatRate,
                                                              const double maximumAllowedHeatLoad,
                                                              const double minimumAllowedDynamicPressure,
                                                              const double minimumAllowedLifetime )
    {
        // Create guidance system object
        guidanceSystem_ = boost::make_shared< GuidanceSystem >( targetPeriapsisAltitude, targetApoapsisAltitude, maximumAllowedHeatRate,
                                                                maximumAllowedHeatLoad, minimumAllowedDynamicPressure, minimumAllowedLifetime,
                                                                true );

        // Finalize guidance system creation
        guidanceSystem_->setCurrentOrbitCounter( 0 );
        guidanceSystem_->createGuidanceSystemObjects(
                    boost::bind( &GuidanceSystemGenerator::propagateTranslationalStateWithCustomTerminationSettings, this, _1, _2, -1.0 ),
                    onboardBodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( ),
                    onboardBodyMap_.at( planetName_ )->getShapeModel( )->getAverageRadius( ) );

        // Give output
        return guidanceSystem_;
    }

    //! Function to propagate state.
    std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > >
    propagateSpacecraftState( const double initialTime, const double finalTime, const Eigen::Vector6d& initialConditions )
    {
        // Set initial time
        integratorSettings_->initialTime_ = initialTime;

        // Set initial state
        propagatorSettings_->resetInitialStates( initialConditions );

        // Set propagation settings
        propagatorSettings_->resetTerminationSettings(
                    boost::make_shared< propagators::PropagationTimeTerminationSettings >( finalTime ) );

        // Create dynamics simulator and propagate
        bool wasPropagationSuccessful = true;
        dynamicsSimulator_ = boost::make_shared< propagators::SingleArcDynamicsSimulator< > >(
                    bodyMap_, integratorSettings_, propagatorSettings_ );
        if ( dynamicsSimulator_->getPropagationTerminationReason( )->getPropagationTerminationReason( ) !=
             propagators::termination_condition_reached )
        {
            wasPropagationSuccessful = false;
        }

        // Extract and give out results
        return std::make_pair( wasPropagationSuccessful, std::make_pair( dynamicsSimulator_->getEquationsOfMotionNumericalSolution( ),
                                                                         dynamicsSimulator_->getDependentVariableHistory( ) ) );
    }

    //! Get current time.
    double getCurrentTime( ) { return currentTime_; }

    //! Get standard gravitational parameter.
    double getGravitationalParameter( )
    {
        return bodyMap_.at( planetName_ )->getGravityFieldModel( )->getGravitationalParameter( );
    }

    //! Get radius.
    double getRadius( )
    {
        return bodyMap_.at( planetName_ )->getShapeModel( )->getAverageRadius( );
    }

    //! Set current time.
    void setCurrentTime( const double currentTime ) { currentTime_ = currentTime; }

private:

    //! Function to propagate translational Cartesian state to specified termination settings.
    std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > >
    propagateTranslationalStateWithCustomTerminationSettings(
            const boost::shared_ptr< propagators::PropagationTerminationSettings > propagationTerminationSettings,
            const Eigen::Vector6d& initialTranslationalCartesianState,
            const double initialTime = -1.0 )
    {
        // Set initial time
        if ( static_cast< int >( initialTime ) == -1 )
        {
            onboardIntegratorSettings_->initialTime_ = currentTime_;
        }
        else
        {
            onboardIntegratorSettings_->initialTime_ = initialTime;
        }

        // Set initial state
        onboardPropagatorSettings_->resetInitialStates( initialTranslationalCartesianState );

        // Set propagation settings
        onboardPropagatorSettings_->resetTerminationSettings( propagationTerminationSettings );

        // Create dynamics simulator and propagate
        bool wasPropagationSuccessful = true;
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator(
                    onboardBodyMap_, onboardIntegratorSettings_, onboardPropagatorSettings_ );
        if ( dynamicsSimulator.getPropagationTerminationReason( )->getPropagationTerminationReason( ) !=
             propagators::termination_condition_reached )
        {
            wasPropagationSuccessful = false;
        }

        // Retrieve results from onboard computer and systems
        std::map< double, Eigen::VectorXd > translationalStateResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariablesResult = dynamicsSimulator.getDependentVariableHistory( );

        // Give output
        return std::make_pair( wasPropagationSuccessful, std::make_pair( translationalStateResult, dependentVariablesResult ) );
    }

    //! Simulation start epoch.
    double simulationStartEpoch_;

    //! Current time.
    double currentTime_;

    //! Name of spacecraft.
    std::string spacecraftName_;

    //! Name of planet.
    std::string planetName_;

    //! Body map pointer.
    tudat::simulation_setup::NamedBodyMap bodyMap_;

    //! Acceleration model map pointer.
    tudat::basic_astrodynamics::AccelerationMap accelerationModelMap_;

    //! Body map pointer.
    tudat::simulation_setup::NamedBodyMap onboardBodyMap_;

    //! Acceleration model map pointer.
    tudat::basic_astrodynamics::AccelerationMap onboardAccelerationModelMap_;

    //! Dynamics simulator pointer.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > >  integratorSettings_;

    //! Dynamics simulator pointer.
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > propagatorSettings_;

    //! Pointer to integrator settings for onboard propagation.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > onboardIntegratorSettings_;

    //! Pointer to propagator settings for onboard propagation.
    boost::shared_ptr< propagators::TranslationalStatePropagatorSettings< > > onboardPropagatorSettings_;

    //! Dynamics simulator pointer.
    boost::shared_ptr< propagators::SingleArcDynamicsSimulator< > > dynamicsSimulator_;

    //! Navigation system pointer.
    boost::shared_ptr< GuidanceSystem > guidanceSystem_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_GENERATORS_H
