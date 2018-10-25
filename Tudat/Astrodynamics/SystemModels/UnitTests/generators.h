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
        const double simulationStartEpoch = 1.0e3;
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
                    boost::make_shared< simulation_setup::SphericalHarmonicAccelerationSettings >( 4, 4 ) );
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

        ///////////////////////////////////////////////  Create instrument model   ///////////////////////////////////////////////////

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
        simulationStartEpoch_ = 1.0e3;
        const double simulationEndEpoch = simulationStartEpoch_ + 1.0;

        // Set body names
        spacecraftName_ = "Satellite";
        planetName_ = "Mars";
        std::vector< std::string > bodiesToCreate = { planetName_ };

        ///////////////////////////////////////////////  Create body map   ///////////////////////////////////////////////////////////////

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

        ///////////////////////////////////////////////  Create acceleration map   ///////////////////////////////////////////////////////

        // Define acceleration settings for simulation model
        std::map< std::string, std::vector< boost::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::SphericalHarmonicAccelerationSettings >( 4, 4 ) );
        accelerationsOfSatellite[ planetName_ ].push_back(
                    boost::make_shared< simulation_setup::AccelerationSettings >( basic_astrodynamics::aerodynamic ) );

        // Set accelerations settings
        simulation_setup::SelectedAccelerationMap accelerationMap;
        accelerationMap[ spacecraftName_ ] = accelerationsOfSatellite;
        onboardAccelerationModelMap_ = simulation_setup::createAccelerationModelsMap(
                    onboardBodyMap_, accelerationMap, { spacecraftName_ }, { planetName_ } );
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
                    frequencyOfDeepSpaceNetworkTracking );

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
        const unsigned int testCase )
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

    // Set simulation time settings.
    const double simulationStartEpoch = 7.0 * tudat::physical_constants::JULIAN_YEAR +
            30.0 * 6.0 * tudat::physical_constants::JULIAN_DAY;
    double simulationEndEpoch = simulationStartEpoch;
    switch ( testCase )
    {
    case 0:
        simulationEndEpoch += 1.5 * tudat::physical_constants::JULIAN_DAY;
        break;
    case 1:
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
    aerodynamicForceCoefficientFiles[ 2 ] = getTudatRootPath( ) + "External/MROLiftCoefficients.txt";

    // Create aerodynamic coefficient settings
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                aerodynamicForceCoefficientFiles, referenceAreaAerodynamic,
                std::vector< AerodynamicCoefficientsIndependentVariables >{ angle_of_attack_dependent, altitude_dependent } );

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
    Eigen::Vector6d initialStateInKeplerianElements;
    switch ( testCase )
    {
    case 0:
    {
        initialStateInKeplerianElements( semiMajorAxisIndex ) = 25946932.3;
        initialStateInKeplerianElements( eccentricityIndex ) = 0.8651912;
        break;
    }
    case 1:
    {
        initialStateInKeplerianElements( semiMajorAxisIndex ) = 4699198.5;
        initialStateInKeplerianElements( eccentricityIndex ) = 0.2546816;
        break;
    }
    default:
        throw std::runtime_error( "Error. Only cases 0 and 1 are supported." );
    }
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 93.0 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 158.7 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 180.0 );

    // Convert to Cartesian coordinates
    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d translationalInitialState = convertKeplerianToCartesianElements(
                initialStateInKeplerianElements, marsGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Dependent variables
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::aerodynamic, "Satellite", "Mars", true ) );
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

} // namespace system_models

} // namespace tudat

#endif // TUDAT_GENERATORS_H
