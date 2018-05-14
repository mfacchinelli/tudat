/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATETORQUEPARTIALS_H
#define TUDAT_CREATETORQUEPARTIALS_H

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EstimationSetup/createAccelerationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/secondDegreeGravitationalTorquePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/sphericalHarmonicGravitationalTorquePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/RotationalDynamicsPartials/inertialTorquePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< acceleration_partials::TorquePartial > createConstantTorqueRotationalDynamicsPartial(
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const basic_astrodynamics::SingleBodyTorqueModelMap& torqueVector );

//! Function to create a single torque partial derivative object.
/*!
 *  Function to create a single torque partial derivative object.
 *  \param torqueModel Torque model for which a partial derivative is to be computed.
 *  \param acceleratedBody Pair of name and object of body undergoing torque
 *  \param acceleratingBody Pair of name and object of body exerting torque
 *  \param bodyMap List of all body objects
 *  \param parametersToEstimate List of parameters that are to be estimated. Empty by default, only required for selected
 *  types of partials (e.g. spherical harmonic torque w.r.t. rotational parameters).
 *  \return Single torque partial derivative object.
 */
template< typename InitialStateParameterType = double >
boost::shared_ptr< acceleration_partials::TorquePartial > createAnalyticalTorquePartial(
        boost::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const basic_astrodynamics::SingleBodyTorqueModelMap& torqueVector = basic_astrodynamics::SingleBodyTorqueModelMap( ),
        const simulation_setup::NamedBodyMap& bodyMap = simulation_setup::NamedBodyMap( ),
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate =
        boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >( ) )
{
    using namespace gravitation;
    using namespace basic_astrodynamics;
    using namespace electro_magnetism;
    using namespace aerodynamics;
    using namespace acceleration_partials;

    boost::shared_ptr< acceleration_partials::TorquePartial > torquePartial;

    // Identify current torque model type
    AvailableTorque torqueType = getTorqueModelType( torqueModel );
    switch( torqueType )
    {
    case second_order_gravitational_torque:

        // Check if identifier is consistent with type.
        if( boost::dynamic_pointer_cast< SecondDegreeGravitationalTorqueModel >( torqueModel ) == NULL )
        {
            throw std::runtime_error( "Torque class type does not match torque type (second_order_gravitational_torque) when making torque partial" );
        }
        else
        {
            // Create partial-calculating object.
            boost::function< double( ) > inertiaTensorNormalizationFunction;
            if( boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >( acceleratedBody.second->getGravityFieldModel( ) )
                    != NULL )
            {
                inertiaTensorNormalizationFunction =
                        boost::bind( &SphericalHarmonicsGravityField::getInertiaTensorNormalizationFactor,
                                     boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                                         acceleratedBody.second->getGravityFieldModel( ) ) );
            }
            torquePartial = boost::make_shared< SecondDegreeGravitationalTorquePartial >
                    ( boost::dynamic_pointer_cast< SecondDegreeGravitationalTorqueModel >( torqueModel ),
                      inertiaTensorNormalizationFunction, acceleratedBody.first, acceleratingBody.first );
        }
        break;
    case spherical_harmonic_gravitational_torque:

        // Check if identifier is consistent with type.
        if( boost::dynamic_pointer_cast< SphericalHarmonicGravitationalTorqueModel >( torqueModel ) == NULL )
        {
            throw std::runtime_error( "Torque class type does not match torque type (spherical_harmonic_gravitational_torque) when making torque partial" );
        }
        else
        {
            // Create partial-calculating object.
            torquePartial = boost::make_shared< SphericalHarmonicGravitationalTorquePartial >
                    ( boost::dynamic_pointer_cast< SphericalHarmonicGravitationalTorqueModel >( torqueModel ),
                      boost::dynamic_pointer_cast< SphericalHarmonicsGravityPartial >( createAnalyticalAccelerationPartial(
                          boost::dynamic_pointer_cast< SphericalHarmonicGravitationalTorqueModel >(
                              torqueModel )->getSphericalHarmonicAcceleration( ), acceleratedBody, acceleratingBody,
                          bodyMap, parametersToEstimate ) ),
                      acceleratedBody.first, acceleratingBody.first );
        }
        break;
    case inertial_torque:
    {
        boost::function< Eigen::Vector3d( ) > angularVelocityFunction =
                boost::bind( &Body::getCurrentAngularVelocityVectorInLocalFrame, acceleratedBody.second );
        boost::function< Eigen::Matrix3d( ) > inertiaTensorFunction =
                boost::bind( &Body::getBodyInertiaTensor, acceleratedBody.second );

        boost::function< double( ) > gravitationalParameterFunction;
        if( acceleratedBody.second->getGravityFieldModel( ) != NULL )
        {
            gravitationalParameterFunction =
                    boost::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                 acceleratedBody.second->getGravityFieldModel( ) );
        }

        boost::function< double( ) > inertiaTensorNormalizationFunction;
        if( boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( acceleratedBody.second->getGravityFieldModel( ) )
                != NULL )
        {
            inertiaTensorNormalizationFunction =
                    boost::bind( &gravitation::SphericalHarmonicsGravityField::getInertiaTensorNormalizationFactor,
                                 boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                     acceleratedBody.second->getGravityFieldModel( ) ) );
        }

        torquePartial = boost::make_shared< acceleration_partials::InertialTorquePartial >(
                    angularVelocityFunction, inertiaTensorFunction, inertiaTensorNormalizationFunction,
                    gravitationalParameterFunction, acceleratedBody.first );
        break;
    }
    default:
        std::string errorMessage = "Torque model " + std::to_string( torqueType ) +
                " not found when making torque partial";
        throw std::runtime_error( errorMessage );
        break;
    }

    return torquePartial;
}

//! This function creates torque partial objects for translational dynamics
/*!
 *  This function creates torque partial objects for translational dynamics from torque models and
 *  list of bodies' states of which derivatives are needed. The return type is an StateDerivativePartialsMap,
 *  a standardized type for communicating such lists of these objects.
 *  \param torqueMap Map of maps containing list of torque models, identifying which torque acts on which
 *   body.
 *  \param bodyMap List of body objects constituting environment for calculations.
 *  \param parametersToEstimate List of parameters which are to be estimated.
 *  \return List of torque-partial-calculating objects in StateDerivativePartialsMap type.
 */
template< typename InitialStateParameterType >
orbit_determination::StateDerivativePartialsMap createTorquePartialsMap(
        const basic_astrodynamics::TorqueModelMap& torqueMap,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate )
{
    // Declare return map.
    orbit_determination::StateDerivativePartialsMap torquePartialsList;

    std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatable_parameters::getListOfRotationalStateParametersToEstimate( parametersToEstimate );
    torquePartialsList.resize( initialDynamicalParameters.size( ) );

    // Iterate over list of bodies of which the partials of the torques acting on them are required.
    for( basic_astrodynamics::TorqueModelMap::const_iterator torqueIterator = torqueMap.begin( );
         torqueIterator != torqueMap.end( ); torqueIterator++ )
    {
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( initialDynamicalParameters.at( i )->getParameterName( ).second.first == torqueIterator->first )
            {
                if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_rotational_body_state ) )
                {
                    // Get object for body undergoing torque
                    const std::string acceleratedBody = torqueIterator->first;
                    boost::shared_ptr< simulation_setup::Body > acceleratedBodyObject = bodyMap.at( acceleratedBody );

                    // Retrieve list of torques acting on current body.
                    basic_astrodynamics::SingleBodyTorqueModelMap torqueVector =
                            torqueMap.at( acceleratedBody );

                    // Declare list of torque partials of current body.
                    std::vector< boost::shared_ptr< orbit_determination::StateDerivativePartial > > torquePartialVector;
                    torquePartialVector.push_back( createConstantTorqueRotationalDynamicsPartial(
                                                      std::make_pair( acceleratedBody, acceleratedBodyObject ),
                                                       torqueVector ) );

                    // Iterate over all torque models and generate their partial-calculating objects.
                    for(  basic_astrodynamics::SingleBodyTorqueModelMap::iterator
                          innerTorqueIterator = torqueVector.begin( );
                          innerTorqueIterator != torqueVector.end( ); innerTorqueIterator++ )
                    {
                        // Get object for body exerting torque
                        std::string acceleratingBody = innerTorqueIterator->first;
                        boost::shared_ptr< simulation_setup::Body > acceleratingBodyObject;
                        if( acceleratingBody != "" )
                        {
                            acceleratingBodyObject = bodyMap.at( acceleratingBody );
                        }

                        for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                        {
                            // Create single partial object
                            boost::shared_ptr< acceleration_partials::TorquePartial > currentTorquePartial =
                                    createAnalyticalTorquePartial(
                                        innerTorqueIterator->second[ j ],
                                        std::make_pair( acceleratedBody, acceleratedBodyObject ),
                                        std::make_pair( acceleratingBody, acceleratingBodyObject ),
                                        torqueVector, bodyMap, parametersToEstimate );

                            torquePartialVector.push_back( currentTorquePartial );
                        }
                    }

                    // Add partials of current body's torques to vector.
                    torquePartialsList[ i ] = torquePartialVector;
                }
            }
        }
    }
    return torquePartialsList;
}

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATETORQUEPARTIALS_H
