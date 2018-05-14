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

#include <limits>
#include <string>
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/numericalAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/gravitationalParameter.h"
#include "Tudat/SimulationSetup/EstimationSetup/createTorquePartials.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electro_magnetism;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_torque_partials )

BOOST_AUTO_TEST_CASE( testSecondDegreeGravitationalTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< ephemerides::ConstantEphemeris >(
                                         boost::lambda::constant( Eigen::Vector6d::Zero( ) ) ) );
    bodyMap[ "Mars" ]->setGravityFieldModel(
                boost::make_shared< gravitation::GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodyMap[ "Phobos" ] = boost::make_shared< Body >( );

    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;


    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );
    bodyMap[ "Phobos" ]->setBodyInertiaTensor(
                phobosInertiaTensor, ( 0.3615 + 0.4265 + 0.5024 ) / 3.0 );

    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
            phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                boost::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed",
                    boost::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true ) ) );

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodyMap[ "Phobos" ]->setRotationalEphemeris( boost::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );



    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;

    bodyMap[ "Phobos" ]->setEphemeris( boost::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                           "Mars", "ECLIPJ2000" ) );


    // Create empty bodies, phobos and mars.
    boost::shared_ptr< Body > phobos = bodyMap.at( "Phobos" );
    boost::shared_ptr< Body > mars = bodyMap.at( "Mars" );
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );

    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
    phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

    mars->setStateFromEphemeris( testTime );
    //    mars->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );


    // Create acceleration due to mars on phobos.
    boost::shared_ptr< SecondDegreeGravitationalTorqueModel > gravitationalTorque =
            createSecondDegreeGravitationalTorqueModel( bodyMap.at( "Phobos" ), bodyMap.at( "Mars" ), "Phobos", "Mars" );


    // Create central gravity partial.
    boost::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( gravitationalTorque, std::make_pair( "Phobos", phobos ),
                                           std::make_pair( "Mars", mars ) );


    // Create gravitational parameter object.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter) );

    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    boost::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodyMap );

    boost::shared_ptr< EstimatableParameter< double > > marsGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 0 );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 1 );
    //    boost::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter = boost::make_shared<
    //            GravitationalParameter >( phobosGravityFieldModel, "Phobos" );

    // Calculate analytical partials.
    torquePartial->update( testTime );

    Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4  );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 4, 3 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtMarsGravitationalParameter = torquePartial->wrtParameter(
                marsGravitationalParameterParameter );
    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter(
                phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter(
                phobosSineCoefficientsParameter );

    //    Eigen::Vector3d partialWrtPhobosGravitationalParameter = centralGravitationPartial->wrtParameter(
    //                phobosGravitationalParameterParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix< double, 3, 4 > testPartialWrtMarsOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );

    Eigen::Matrix3d testPartialWrtMarsPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMarsVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in orientation for numerical partial/
    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation <<  1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    boost::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, _1 );
    boost::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, _1 );

    //    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange(
                phobosRotationalStateSetFunction, gravitationalTorque, phobos->getRotationalStateVector( ), orientationPerturbation,
                appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, gravitationalTorque, phobos->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    testPartialWrtMarsOrientation =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, gravitationalTorque, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, gravitationalTorque, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );



    boost::function< void( Eigen::Vector6d ) > phobosStateSetFunction =
            boost::bind( &Body::setState,  phobos, _1 );
    boost::function< void( Eigen::Vector6d ) > marsStateSetFunction =
            boost::bind( &Body::setState, mars, _1 );

    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 100.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Calculate numerical partials.
    testPartialWrtMarsPosition = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, gravitationalTorque, mars->getState( ), positionPerturbation, 0 );
    testPartialWrtMarsVelocity = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, gravitationalTorque, mars->getState( ), velocityPerturbation, 3 );
    testPartialWrtPhobosPosition = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, gravitationalTorque, phobos->getState( ), positionPerturbation, 0 );
    testPartialWrtPhobosVelocity = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, gravitationalTorque, phobos->getState( ), velocityPerturbation, 3 );

    Eigen::Vector3d testPartialWrtMarsGravitationalParameter = calculateTorqueWrtParameterPartials(
                marsGravitationalParameterParameter, gravitationalTorque, 1.0E12 );

    boost::function< void( ) > updateFunction = &emptyFunction;
            //boost::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
                phobosCosineCoefficientsParameter, gravitationalTorque,
                Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateTorqueWrtParameterPartials(
                phobosSineCoefficientsParameter, gravitationalTorque,
                Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque,
                                           analyticalChangeInTorque, 1.0E-6 );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity,
                                       partialWrtPhobosRotationalVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsOrientation,
                                       partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationalVelocity,
                                       partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );


    //    std::cout<<partialWrtPhobosState<<std::endl<<std::endl<<
    //               partialWrtMarsState<<std::endl<<std::endl;
    //    std::cout<<testPartialWrtMarsPosition<<std::endl<<std::endl<<
    //               testPartialWrtMarsVelocity<<std::endl<<std::endl<<
    //               testPartialWrtPhobosPosition<<std::endl<<std::endl<<
    //               testPartialWrtPhobosVelocity<<std::endl<<std::endl;

    //    std::cout<<( partialWrtPhobosState.block( 0, 0, 3, 3 ) - testPartialWrtPhobosPosition ).cwiseQuotient(
    //                   testPartialWrtPhobosPosition )<<std::endl<<std::endl;
    //    std::cout<<( partialWrtMarsState.block( 0, 0, 3, 3 ) - testPartialWrtMarsPosition ).cwiseQuotient(
    //                   testPartialWrtMarsPosition )<<std::endl<<std::endl;


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 0, 3, 3 ),
                                        testPartialWrtPhobosPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 0, 3, 3 ),
                                        testPartialWrtMarsPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 3, 3, 3 ),
                                        testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 3, 3, 3 ),
                                        testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsGravitationalParameter,
                                       partialWrtMarsGravitationalParameter, 1.0E-6 );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 2, 2 ) ), 1.0E6 );
    testPartialWrtPhobosCosineCoefficients( 2, 2 ) = 0.0;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosCosineCoefficients,
                                       partialWrtPhobosCosineCoefficients, 1.0E-9 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosSineCoefficients,
                                       partialWrtPhobosSineCoefficients, 1.0E-9 );


}

BOOST_AUTO_TEST_CASE( testInertialTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< ephemerides::ConstantEphemeris >(
                                         boost::lambda::constant( Eigen::Vector6d::Zero( ) ) ) );
    bodyMap[ "Mars" ]->setGravityFieldModel(
                boost::make_shared< gravitation::GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodyMap[ "Phobos" ] = boost::make_shared< Body >( );

    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;


    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );
    bodyMap[ "Phobos" ]->setBodyInertiaTensor(
                phobosInertiaTensor, ( 0.3615 + 0.4265 + 0.5024 ) / 3.0 );

    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
            phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                boost::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed" ) );

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );
    unitRotationState( 4 ) = 1.0E-5;
    unitRotationState( 5 ) = 2.0E-5;
    unitRotationState( 6 ) = -3.2E-5;

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodyMap[ "Phobos" ]->setRotationalEphemeris( boost::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );



    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;

    bodyMap[ "Phobos" ]->setEphemeris( boost::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                           "Mars", "ECLIPJ2000" ) );


    // Create empty bodies, phobos and mars.
    boost::shared_ptr< Body > phobos = bodyMap.at( "Phobos" );
    boost::shared_ptr< Body > mars = bodyMap.at( "Mars" );
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );

    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
    phobosRotationalState( 4 ) = 1.0E-5;
    phobosRotationalState( 5 ) = 2.0E-5;
    phobosRotationalState( 6 ) = -3.2E-5;

    phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

    mars->setStateFromEphemeris( testTime );
    //    mars->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );


    // Create acceleration due to mars on phobos.
    boost::shared_ptr< InertialTorqueModel > inertialTorqueModel =
            createInertialTorqueModel( bodyMap.at( "Phobos" ), "Phobos" );
    inertialTorqueModel->updateMembers( 0.0 );

    SingleBodyTorqueModelMap torqueList;
    torqueList[ "Phobos" ].push_back( inertialTorqueModel );

    // Create central gravity partial.
    boost::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( inertialTorqueModel, std::make_pair( "Phobos", phobos ),
                                           std::make_pair( "Phobos", phobos ) );


    // Create gravitational parameter object.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );

    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    boost::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodyMap );

    boost::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    boost::shared_ptr< EstimatableParameter< double > > meanMomentOfInertiaParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 1 );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 0 );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 1 );
    //    boost::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter = boost::make_shared<
    //            GravitationalParameter >( phobosGravityFieldModel, "Phobos" );

    // Calculate analytical partials.
    torquePartial->update( testTime );

    Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4  );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 4, 3 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtPhobosGravitationalParameter = torquePartial->wrtParameter(
                phobosGravitationalParameterParameter );
    Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter(
                meanMomentOfInertiaParameter );

    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter(
                phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter(
                phobosSineCoefficientsParameter );


    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix< double, 3, 4 > testPartialWrtMarsOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );

    Eigen::Matrix3d testPartialWrtMarsPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMarsVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in orientation for numerical partial/
    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation <<  1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    boost::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, _1 );
    boost::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, _1 );

    //    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange(
                phobosRotationalStateSetFunction, inertialTorqueModel, phobos->getRotationalStateVector( ), orientationPerturbation,
                appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, inertialTorqueModel, phobos->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    testPartialWrtMarsOrientation =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, inertialTorqueModel, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, inertialTorqueModel, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );



    boost::function< void( Eigen::Vector6d ) > phobosStateSetFunction =
            boost::bind( &Body::setState,  phobos, _1 );
    boost::function< void( Eigen::Vector6d ) > marsStateSetFunction =
            boost::bind( &Body::setState, mars, _1 );

    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 100.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Calculate numerical partials.
    testPartialWrtMarsPosition = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, inertialTorqueModel, mars->getState( ), positionPerturbation, 0 );
    testPartialWrtMarsVelocity = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, inertialTorqueModel, mars->getState( ), velocityPerturbation, 3 );
    testPartialWrtPhobosPosition = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, inertialTorqueModel, phobos->getState( ), positionPerturbation, 0 );
    testPartialWrtPhobosVelocity = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, inertialTorqueModel, phobos->getState( ), velocityPerturbation, 3 );


    boost::function< void( ) > updateFunction = //&emptyFunction;
            boost::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true );
    Eigen::Vector3d testPartialWrtPhobosGravitationalParameter = calculateTorqueWrtParameterPartials(
                phobosGravitationalParameterParameter, inertialTorqueModel, 1.0E8, updateFunction );
    Eigen::MatrixXd testPartialWrtMeanMomentOfInertia = calculateTorqueWrtParameterPartials(
                meanMomentOfInertiaParameter, inertialTorqueModel, 1.0E-4, updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
                phobosCosineCoefficientsParameter, inertialTorqueModel,
                Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateTorqueWrtParameterPartials(
                phobosSineCoefficientsParameter, inertialTorqueModel,
                Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque,
                                           analyticalChangeInTorque, std::numeric_limits< double >::epsilon( ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity,
                                       partialWrtPhobosRotationalVelocity, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsOrientation,
                                       partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationalVelocity,
                                       partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 0, 3, 3 ),
                                        testPartialWrtPhobosPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 0, 3, 3 ),
                                        testPartialWrtMarsPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 3, 3, 3 ),
                                        testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 3, 3, 3 ),
                                        testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosGravitationalParameter,
                                       partialWrtPhobosGravitationalParameter, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 0, 3, 2 ),
                                       Eigen::MatrixXd::Zero( 3, 2 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 5, 3, 4 ),
                                       Eigen::MatrixXd::Zero( 3, 4 ), std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 0, 3, 1 ),
                                       Eigen::MatrixXd::Zero( 3, 1 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 3, 3, 3 ),
                                       Eigen::MatrixXd::Zero( 3, 3 ), std::numeric_limits< double >::epsilon( ) );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 2, 2 ) ), 1.0E5 );
    testPartialWrtPhobosCosineCoefficients( 2, 2 ) = 0.0;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosCosineCoefficients,
                                       partialWrtPhobosCosineCoefficients, 1.0E-9 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosSineCoefficients,
                                       partialWrtPhobosSineCoefficients, 1.0E-9 );

    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( partialWrtMeanMomentOfInertia( i, 0 ) ), 1.0E3 );
    }
}

class EffectiveTorqueModel: public basic_astrodynamics::TorqueModel
{
public:
    EffectiveTorqueModel(
            const boost::function< Eigen::Matrix3d( ) > inertiaTensorFunction,
            const SingleBodyTorqueModelMap& torqueList ):
        inertiaTensorFunction_( inertiaTensorFunction ), torqueList_( torqueList ){ }


    Eigen::Vector3d getTorque( )
    {
        return effectiveTorque_;
    }

    void updateMembers( const double currentTime )
    {
        Eigen::Vector3d totalTorque = Eigen::Vector3d::Zero( );
        for( auto it = torqueList_.begin( ); it != torqueList_.end( ); it++ )
        {
            for( unsigned int i = 0; i < it->second.size( ); i++ )
            {
                totalTorque += it->second.at( i )->getTorque( );
            }
        }

        effectiveTorque_ = inertiaTensorFunction_( ).inverse( ) * totalTorque;
    }

private:

    boost::function< Eigen::Matrix3d( ) > inertiaTensorFunction_;

    SingleBodyTorqueModelMap torqueList_;

    Eigen::Vector3d effectiveTorque_;
};

BOOST_AUTO_TEST_CASE( testConstantTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< ephemerides::ConstantEphemeris >(
                                         boost::lambda::constant( Eigen::Vector6d::Zero( ) ) ) );
    bodyMap[ "Mars" ]->setGravityFieldModel(
                boost::make_shared< gravitation::GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodyMap[ "Phobos" ] = boost::make_shared< Body >( );

    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;


    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );
    bodyMap[ "Phobos" ]->setBodyInertiaTensor(
                phobosInertiaTensor, ( 0.3615 + 0.4265 + 0.5024 ) / 3.0 );

    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
            phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                boost::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed",
                    boost::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true ) ) );

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );
    unitRotationState( 4 ) = 1.0E-5;
    unitRotationState( 5 ) = 2.0E-5;
    unitRotationState( 6 ) = -3.2E-5;

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodyMap[ "Phobos" ]->setRotationalEphemeris( boost::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );



    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;

    bodyMap[ "Phobos" ]->setEphemeris( boost::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                           "Mars", "ECLIPJ2000" ) );


    // Create empty bodies, phobos and mars.
    boost::shared_ptr< Body > phobos = bodyMap.at( "Phobos" );
    boost::shared_ptr< Body > mars = bodyMap.at( "Mars" );
    setGlobalFrameBodyEphemerides( bodyMap, "Mars", "ECLIPJ2000" );

    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );

    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
    phobosRotationalState( 4 ) = 1.0E-5;
    phobosRotationalState( 5 ) = 2.0E-5;
    phobosRotationalState( 6 ) = -3.2E-5;

    phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

    mars->setStateFromEphemeris( testTime );
    //    mars->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );


    // Create acceleration due to mars on phobos.
    boost::shared_ptr< SecondDegreeGravitationalTorqueModel > gravitationalTorque =
            createSecondDegreeGravitationalTorqueModel( bodyMap.at( "Phobos" ), bodyMap.at( "Mars" ), "Phobos", "Mars" );
    boost::shared_ptr< InertialTorqueModel > inertialTorqueModel =
            createInertialTorqueModel( bodyMap.at( "Phobos" ), "Phobos" );
    gravitationalTorque->updateMembers( 0.0 );
    inertialTorqueModel->updateMembers( 0.0 );

    SingleBodyTorqueModelMap torqueList;
    torqueList[ "Mars" ].push_back( gravitationalTorque );
    torqueList[ "Phobos" ].push_back( inertialTorqueModel );

    boost::shared_ptr< EffectiveTorqueModel > effectiveTorqueModel =
            boost::make_shared< EffectiveTorqueModel >(
                boost::bind( &Body::getBodyInertiaTensor, bodyMap.at( "Phobos" ) ), torqueList );
    effectiveTorqueModel->updateMembers( 0.0 );

    // Create central gravity partial.
    boost::shared_ptr< TorquePartial > torquePartial =
            createConstantTorqueRotationalDynamicsPartial(
                std::make_pair( "Phobos", bodyMap.at( "Phobos" ) ), torqueList );

    // Create gravitational parameter object.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );

    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    boost::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodyMap );

    boost::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    boost::shared_ptr< EstimatableParameter< double > > meanMomentOfInertiaParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 1 );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 0 );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 1 );

    // Calculate analytical partials.
    torquePartial->update( testTime );

    Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4  );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 4, 3 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtPhobosGravitationalParameter = torquePartial->wrtParameter(
                phobosGravitationalParameterParameter );
    Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter(
                meanMomentOfInertiaParameter );

    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter(
                phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter(
                phobosSineCoefficientsParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix< double, 3, 4 > testPartialWrtMarsOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );

    Eigen::Matrix3d testPartialWrtMarsPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMarsVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in orientation for numerical partial/
    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation <<  1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    boost::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, _1 );
    boost::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, _1 );

    //    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange(
                phobosRotationalStateSetFunction, effectiveTorqueModel, phobos->getRotationalStateVector( ), orientationPerturbation,
                appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, effectiveTorqueModel, phobos->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    testPartialWrtMarsOrientation =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, effectiveTorqueModel, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, effectiveTorqueModel, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );



    boost::function< void( Eigen::Vector6d ) > phobosStateSetFunction =
            boost::bind( &Body::setState,  phobos, _1 );
    boost::function< void( Eigen::Vector6d ) > marsStateSetFunction =
            boost::bind( &Body::setState, mars, _1 );

    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 100.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Calculate numerical partials.
    testPartialWrtMarsPosition = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, effectiveTorqueModel, mars->getState( ), positionPerturbation, 0 );
    testPartialWrtMarsVelocity = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, effectiveTorqueModel, mars->getState( ), velocityPerturbation, 3 );
    testPartialWrtPhobosPosition = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, effectiveTorqueModel, phobos->getState( ), positionPerturbation, 0 );
    testPartialWrtPhobosVelocity = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, effectiveTorqueModel, phobos->getState( ), velocityPerturbation, 3 );


    boost::function< void( ) > updateFunction = &emptyFunction;
           // boost::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodyMap.at( "Phobos" ), true );
    Eigen::Vector3d testPartialWrtPhobosGravitationalParameter = calculateTorqueWrtParameterPartials(
                phobosGravitationalParameterParameter, effectiveTorqueModel, 1.0E2, updateFunction );
    Eigen::MatrixXd testPartialWrtMeanMomentOfInertia = calculateTorqueWrtParameterPartials(
                meanMomentOfInertiaParameter, effectiveTorqueModel, 1.0E-4, updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
                phobosCosineCoefficientsParameter, effectiveTorqueModel,
                Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateTorqueWrtParameterPartials(
                phobosSineCoefficientsParameter, effectiveTorqueModel,
                Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque,
                                           analyticalChangeInTorque, std::numeric_limits< double >::epsilon( ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity,
                                       partialWrtPhobosRotationalVelocity, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsOrientation,
                                       partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationalVelocity,
                                       partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 0, 3, 3 ),
                                        testPartialWrtPhobosPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 0, 3, 3 ),
                                        testPartialWrtMarsPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 3, 3, 3 ),
                                        testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 3, 3, 3 ),
                                        testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( Eigen::MatrixXd( phobosInertiaTensor ) * testPartialWrtPhobosGravitationalParameter ),
                                       partialWrtPhobosGravitationalParameter, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( Eigen::MatrixXd( phobosInertiaTensor ) * testPartialWrtMeanMomentOfInertia ),
                                       partialWrtMeanMomentOfInertia, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 0, 3, 2 ),
                                       Eigen::MatrixXd::Zero( 3, 2 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 5, 3, 4 ),
                                       Eigen::MatrixXd::Zero( 3, 4 ), std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 0, 3, 1 ),
                                       Eigen::MatrixXd::Zero( 3, 1 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 3, 3, 3 ),
                                       Eigen::MatrixXd::Zero( 3, 3 ), std::numeric_limits< double >::epsilon( ) );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    BOOST_CHECK_SMALL( std::fabs(
                           ( Eigen::MatrixXd( phobosInertiaTensor ) * testPartialWrtPhobosCosineCoefficients )( 2, 4 ) ), 1.0E5 );
    testPartialWrtPhobosCosineCoefficients( 2, 4 ) = 0.0;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::MatrixXd( phobosInertiaTensor ) * testPartialWrtPhobosCosineCoefficients ), partialWrtPhobosCosineCoefficients, 1.0E-9 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( Eigen::MatrixXd( phobosInertiaTensor ) * testPartialWrtPhobosSineCoefficients ), partialWrtPhobosSineCoefficients, 1.0E-9 );

}

BOOST_AUTO_TEST_SUITE_END( )


} // namespace unit_tests

} // namespace tudat




