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

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/orbitDeterminationManager.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_multi_arc_state_estimation )

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::unit_conversions;

template< typename ObservationScalarType = double , typename TimeType = double , typename StateScalarType  = double >
Eigen::VectorXd  executeParameterEstimation(
        const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference = Eigen::Matrix< StateScalarType, 12, 1 >::Zero( ),
        const Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( 2 ),
        const bool propagateVariationalEquations = 1,
        const bool patchMultiArcs = 0,
        const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > forcedMultiArcInitialStates =
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( ),
        const double arcDuration = 4.0 * 86400.0,
        const double arcOverlap  = 5.0E3 )
{
    //Load spice kernels.f
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Earth" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 20.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Orbiter" ] = boost::make_shared< Body >( );
    bodyMap[ "Orbiter" ]->setConstantBodyMass( 5.0E3 );
    bodyMap[ "Orbiter" ]->setEphemeris( boost::make_shared< MultiArcEphemeris >(
                                            std::map< double, boost::shared_ptr< Ephemeris > >( ),
                                            "Mars", "ECLIPJ2000" ) );

    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > orbiterRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Orbiter" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    orbiterRadiationPressureSettings, "Orbiter", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );

    createGroundStation( bodyMap.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( grazStation );
    groundStations.push_back( mslStation );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap singleArcAccelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    singleArcAccelerationMap[ "Mars" ] = accelerationsOfMars;

    std::vector< std::string > singleArcBodiesToIntegrate, singleArcCentralBodies;
    singleArcBodiesToIntegrate.push_back( "Mars" );
    singleArcCentralBodies.push_back( "SSB" );

    AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, singleArcAccelerationMap, singleArcBodiesToIntegrate, singleArcCentralBodies );
    Eigen::VectorXd singleArcInitialStates = getInitialStatesOfBodies(
                singleArcBodiesToIntegrate, singleArcCentralBodies, bodyMap, initialEphemerisTime );

    singleArcInitialStates += initialStateDifference.segment(
                0, singleArcInitialStates.rows( ) );

    boost::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< > >(
                singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToIntegrate,
                singleArcInitialStates, finalEphemerisTime );


    SelectedAccelerationMap multiArcAccelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfOrbiter;
    accelerationsOfOrbiter[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfOrbiter[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    multiArcAccelerationMap[ "Orbiter" ] = accelerationsOfOrbiter;

    std::vector< std::string > multiArcBodiesToIntegrate, multiArcCentralBodies;
    multiArcBodiesToIntegrate.push_back( "Orbiter" );
    multiArcCentralBodies.push_back( "Mars" );

    AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, multiArcAccelerationMap, multiArcBodiesToIntegrate, multiArcCentralBodies );


    // Creater arc times
    std::vector< double > integrationArcStarts, integrationArcEnds, integrationArcLimits;

    double integrationStartTime = initialEphemerisTime;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double currentStartTime = integrationStartTime;
    double currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcStarts.push_back( currentStartTime );
        integrationArcLimits.push_back( currentStartTime );

        integrationArcEnds.push_back( currentEndTime );

        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );

    integrationArcLimits.push_back( currentStartTime + arcOverlap );

    // Create list of multi-arc initial states
    unsigned int numberOfIntegrationArcs = integrationArcStarts.size( );
    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    multiArcSystemInitialStates.resize( numberOfIntegrationArcs );

    // Define (quasi-arbitrary) arc initial states
    double marsGravitationalParameter =  bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );

    if( forcedMultiArcInitialStates.size( ) == 0 )
    {
        for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
        {
            Eigen::Vector6d orbiterInitialStateInKeplerianElements;
            orbiterInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6000.0E3;
            orbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
            orbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
            orbiterInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                    = convertDegreesToRadians( 235.7 - j );
            orbiterInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                    = convertDegreesToRadians( 23.4 + j );
            orbiterInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 + j * 10.0 );

            // Convert state from Keplerian elements to Cartesian elements.
            multiArcSystemInitialStates[ j ]  = convertKeplerianToCartesianElements(
                        orbiterInitialStateInKeplerianElements,
                        marsGravitationalParameter ) + initialStateDifference.segment(
                        singleArcInitialStates.rows( ), 6 ).template cast< double >( );
        }
    }
    else
    {
        for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
        {
            multiArcSystemInitialStates[ j ] = forcedMultiArcInitialStates.at( j ) + initialStateDifference.segment(
                        singleArcInitialStates.rows( ), 6 );
        }
    }

    // Create propagation settings for each arc
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
    {
        arcPropagationSettingsList.push_back(
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( multiArcCentralBodies, multiArcAccelerationModelMap, multiArcBodiesToIntegrate,
                      multiArcSystemInitialStates.at( i ), integrationArcEnds.at( i ) ) );
    }

    boost::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
            boost::make_shared< MultiArcPropagatorSettings< > >( arcPropagationSettingsList, patchMultiArcs );

    boost::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings =
            boost::make_shared< HybridArcPropagatorSettings< > >(
                singleArcPropagatorSettings, multiArcPropagatorSettings );


    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 60.0 );


    // Set parameters that are to be estimated.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                boost::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                    multiArcBodiesToIntegrate.at( 0 ), multiArcPropagatorSettings->getInitialStates( ),
                    integrationArcStarts, multiArcCentralBodies.at( 0 ) ) );
    parameterNames.push_back(
                boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                    singleArcBodiesToIntegrate.at( 0 ), singleArcInitialStates, singleArcCentralBodies.at( 0 ) ) );

    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );

    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodyMap );


    // Define links in simulation.
    std::vector< LinkEnds > linkEnds2;
    linkEnds2.resize( 4 );
    linkEnds2[ 0 ][ transmitter ] = grazStation;
    linkEnds2[ 0 ][ receiver ] = std::make_pair( "Orbiter", "" );

    linkEnds2[ 1 ][ receiver ] = grazStation;
    linkEnds2[ 1 ][ transmitter ] = std::make_pair( "Orbiter", "" );

    linkEnds2[ 2 ][ transmitter ] = grazStation;
    linkEnds2[ 2 ][ receiver ] = std::make_pair( "Mars", "" );

    linkEnds2[ 3 ][ receiver ] = grazStation;
    linkEnds2[ 3 ][ transmitter ] = std::make_pair( "Mars", "" );

    observation_models::ObservationSettingsMap observationSettingsMap;
    for( unsigned int i = 0; i  < 4; i++ )
    {
        observationSettingsMap.insert( std::make_pair( linkEnds2[ i ], boost::make_shared< ObservationSettings >(
                                                           one_way_range ) ) );
        observationSettingsMap.insert( std::make_pair( linkEnds2[ i ], boost::make_shared< ObservationSettings >(
                                                           angular_position ) ) );
    }


    // Create orbit determination object.
    OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< ObservationScalarType, TimeType >(
                bodyMap, parametersToEstimate,
                observationSettingsMap, integratorSettings, hybridArcPropagatorSettings );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );


    TimeType observationTime;
    int numberOfObservationsPerArc = 5000;
    double timeBuffer = 9000.0;


    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservationsPerArc * integrationArcStarts.size( ) );

    for( unsigned int i = 0; i < integrationArcLimits.size( ) - 1; i++ )
    {
        double currentTimeStep = ( integrationArcLimits[ i + 1 ] - integrationArcLimits[ i ] - 2.0 * timeBuffer ) /
                static_cast< double >( numberOfObservationsPerArc - 1 );
        observationTime = integrationArcLimits[ i ] + timeBuffer;
        for( int j = 0; j < numberOfObservationsPerArc; j++ )
        {
            initialObservationTimes[ j + i * numberOfObservationsPerArc ] = observationTime;
            observationTime += currentTimeStep;
        }
    }


    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > measurementSimulationInput;

    for( unsigned int i = 0; i < 4; i++ )
    {
        measurementSimulationInput[ one_way_range ][ linkEnds2[ i ] ] = std::make_pair( initialObservationTimes, receiver );
        measurementSimulationInput[ angular_position ][ linkEnds2[ i ] ] = std::make_pair( initialObservationTimes, receiver );
    }


    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    PodInputDataType observationsAndTimes = simulateObservations< ObservationScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( )  );
    std::cout<<"Number of observables "<<observationsAndTimes.size( )<<std::endl;
    std::cout<<"Number of range link ends  "<<observationsAndTimes.at( one_way_range ).size( )<<std::endl;


    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    std::cout << "Truth " << std::setprecision( 6 ) << truthParameters.transpose( ) << std::endl;

    initialParameterEstimate[ 0 ] += 1.0E2;
    initialParameterEstimate[ 1 ] += 1.0E2;
    initialParameterEstimate[ 2 ] += 1.0E2;
    initialParameterEstimate[ 3 ] += 1E-3;
    initialParameterEstimate[ 4 ] += 1E-3;
    initialParameterEstimate[ 5 ] += 1E-3;

    for( unsigned int i = 1; i < multiArcBodiesToIntegrate.size( ) * integrationArcStarts.size( ); i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 1E-3;
        initialParameterEstimate[ 4 + 6 * i ] += 1E-3;
        initialParameterEstimate[ 5 + 6 * i ] += 1E-3;
    }

    for( unsigned int i = 6 * ( 1 + multiArcBodiesToIntegrate.size( ) * integrationArcStarts.size( ) );
         i < static_cast< unsigned int >( initialParameterEstimate.rows( ) ); i++ )
    {
        initialParameterEstimate[ i ] *= ( 1.0 + 1.0E-6 );
    }

    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    boost::shared_ptr< PodInput< ObservationScalarType, TimeType > > podInput =
            boost::make_shared< PodInput< ObservationScalarType, TimeType > >(
                observationsAndTimes, ( initialParameterEstimate ).rows( ) );

    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0E-4;
    weightPerObservable[ angular_position ] = 1.0E-20;

    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );

    boost::shared_ptr< PodOutput< StateScalarType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput );

    return ( podOutput->parameterEstimate_ - truthParameters ).template cast< double >( );
}


BOOST_AUTO_TEST_CASE( test_MultiArcStateEstimation )
{
    Eigen::VectorXd parameterError = executeParameterEstimation< double, double, double >( );
    int numberOfEstimatedArcs = ( parameterError.rows( ) - 8 ) / 6;

    std::cout << parameterError.transpose( ) << std::endl;

    for( unsigned int j = 0; j < 2; j++ )
    {
        BOOST_CHECK_SMALL( std::fabs( parameterError( j ) ), 1.0 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( j + 3 ) ), 1.0E-6  );
    }

    BOOST_CHECK_SMALL( std::fabs( parameterError( 2 ) ), 500.0 );
    BOOST_CHECK_SMALL( std::fabs( parameterError( 5 ) ), 1.0E-4  );

    for( int i = 0; i < numberOfEstimatedArcs; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( parameterError( ( i + 1 )* 6 + j ) ), 1E-1 );
            BOOST_CHECK_SMALL( std::fabs( parameterError( ( i + 1 ) * 6 + j + 3 ) ), 1.0E-5  );
        }
    }

    BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 2 ) ), 1.0E11 );
    BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 1 ) ), 1.0E6 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
