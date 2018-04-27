/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NUNIFIEDSTATEMODELMODIFIEDRODRIGUESPARAMETERSSTATEDERIVATIVE_H
#define TUDAT_NUNIFIEDSTATEMODELMODIFIEDRODRIGUESPARAMETERSSTATEDERIVATIVE_H

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluate the equations of motion for the unifies state model with modified rodrigues parameters (USM6)
/*!
 * Function to evaluate the equations of motion for the unifies state model with modified rodrigues parameters (USM6), providing the
 * time-derivatives of USM6 elements from the accelerations expressed in an RSW frame (see Vallado, 2001).
 * This function takes a number of precomputed quantities as input, to reduce computational burden
 * \param currentUnifiedStateModelElements Current USM6 elements of the body for which the equations of motion are
 * to be evaluated
 * \param accelerationsInRswFrame Accelerations acting on body, expressed in RSW frame
 * \param sineLambdaParameter Sine of the right ascension of latitude
 * \param cosineLambdaParameter Cosine of the right ascension of latitude
 * \param gammaParameter Value of the parameter gamma (see Vittaldev, 2010)
 * \param rotationalVelocityVector Rotational velocity of the local orbital frame w.r.t. the inertial frame (see Vittaldev, 2010)
 * \param pParameterVector Value of the vector gamma (see Vittaldev, 2010)
 * \return Time derivatives of USM6 elements.
 */
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelModifiedRodriguesParameters(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double sineLambdaParameter,
        const double cosineLambdaParameter,
        const double gammaParameter,
        const Eigen::Vector3d rotationalVelocityVector,
        const Eigen::Vector3d pParameterVector );

//! Function to evaluate the equations of motion for the unifies state model with modified rodrigues parameters (USM6)
/*!
 * Function to evaluate the equations of motion for the unifies state model with modified rodrigues parameters (USM6), providing the
 * time-derivatives of USM6 elements from the accelerations expressed in an RSW frame (see Vallado, 2001).
 * \param currentUnifiedStateModelElements Current USM6 elements of the body for which the equations of motion are
 * to be evaluated
 * \param accelerationsInRswFrame Accelerations acting on body, expressed in RSW frame
 * \param centralBodyGravitationalParameter Gravitational parameter of sum of central body and body for which orbit is propagated.
 * \return Time derivatives of USM6 elements.
 */
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelModifiedRodriguesParameters(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector3d& accelerationsInRswFrame,
        const double centralBodyGravitationalParameter );

//! Function to evaluate the equations of motion for the unifies state model with modified rodrigues parameters (USM6)
/*!
 * Function to evaluate the equations of motion for the unifies state model with modified rodrigues parameters (USM6), providing the
 * time-derivatives of USM6 elements from the accelerations expressed in an RSW frame (see Vallado, 2001). This function takes the accelerations
 * in the inertial frame, as well as the Cartesian inertial state, and converts the accelerations to the RSW frame.
 * \param currentUnifiedStateModelElements Current USM6 elements of the body for which the equations of motion are
 * to be evaluated
 * \param currentCartesianState Current Cartesian state of the body for which the equations of motion are to be evaluated
 * \param accelerationsInInertialFrame Accelerations acting on body, expressed in inertial frame
 * \param centralBodyGravitationalParameter Gravitational parameter of sum of central body and body for which orbit is propagated.
 * \return Time derivatives of USM6 elements.
 */
Eigen::Vector7d computeStateDerivativeForUnifiedStateModelModifiedRodriguesParameters(
        const Eigen::Vector7d& currentUnifiedStateModelElements,
        const Eigen::Vector6d& currentCartesianState,
        const Eigen::Vector3d& accelerationsInInertialFrame,
        const double centralBodyGravitationalParameter );

//! Class for computing the state derivative of translational motion of N bodies, using a Gauss method with
//! unified state model with modified rodrigues parameters (USM6).
/*!
 * Class for computing the state derivative of translational motion of N bodies, using a Gauss method with unified
 * state model with modified rodrigues parameters (USM6). In this method, the derivative of the USM6 elements are computed from the total
 * Cartesian accelerations, with the USM6 elements of the bodies the states being numerically propagated.
 */
template< typename StateScalarType = double, typename TimeType = double >
class NBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on. This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative(
            const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
            const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
            const std::vector< std::string >& bodiesToIntegrate ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, unified_state_model_modified_rodrigues_parameters, bodiesToIntegrate )
    {
        originalAccelerationModelsPerBody_ = this->accelerationModelsPerBody_ ;

        // Remove central gravitational acceleration from list of accelerations that is to be evaluated
        centralBodyGravitationalParameters_ =
                removeCentralGravityAccelerations(
                    centralBodyData->getCentralBodies( ), this->bodiesToBeIntegratedNumerically_,
                    this->accelerationModelsPerBody_ );
        this->createAccelerationModelList( );

    }

    //! Destructor
    ~NBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system, using the equations of motion for the
    //! unified state model with modified rodrigues parameters (USM6).
    /*!
     *  Calculates the state derivative of the translational motion of the system, using the equations of motion for the
     *  unified state model with modified rodrigues parameters (USM6). The input is the current state in USM6 elememts. The state derivate
     *  of this set is computed. To do so the accelerations are internally transformed into the RSW frame, using the current
     *  Cartesian state as set by the last call to the convertToOutputSolution function
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 7 * bodiesToBeIntegratedNumerically_.size( ), containing USM6
     *  elements of the bodies being integrated.
     *  The order of the values is defined by the order of bodies in bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current derivative of the USM6 elements of the
     *  system of bodies integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        // Get total inertial accelerations acting on bodies
        stateDerivative.setZero( );
        this->sumStateDerivativeContributions( stateOfSystemToBeIntegrated, stateDerivative, false );

        // Compute RSW accelerations for each body, and evaluate equations of motion for USM6 elements.
        Eigen::Vector3d currentAccelerationInRswFrame;
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentAccelerationInRswFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                        currentCartesianLocalSoluton_.segment( i * 6, 6 ).template cast< double >( ) ) *
                    stateDerivative.block( i * 6 + 3, 0, 6, 1 ).template cast< double >( );

            stateDerivative.block( i * 7, 0, 7, 1 ) = computeStateDerivativeForUnifiedStateModelModifiedRodriguesParameters(
                        stateOfSystemToBeIntegrated.block( i * 7, 0, 7, 1 ).template cast< double >( ), currentAccelerationInRswFrame,
                        centralBodyGravitationalParameters_.at( i )( ) ).template cast< StateScalarType >( );
        }

    }

    //! Function to convert the state in the conventional form to the USM6 elements form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the USM6 propagator,
     * this transforms the Cartesian state w.r.t. the central body (conventional form) to the UMS7 elements
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid, used to computed Kepler orbits
     * \return State (outputSolution), converted to the USM6 elements
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        // Subtract frame origin and Keplerian states from inertial state.
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentState =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                    this->bodiesToBeIntegratedNumerically_.size( ) * 7 );

        // Convert state to USM6 for each body
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentState.segment( i * 7, 7 ) =
                    orbital_element_conversions::convertCartesianToUnifiedStateModelModifiedRodriguesParametersElements(
                        cartesianSolution.block( i * 6, 0, 6, 1 ).template cast< double >( ), static_cast< double >(
                            centralBodyGravitationalParameters_.at( i )( ) ) ).template cast< StateScalarType >( );
        }

        return currentState;
    }

    //! Function to convert the USM6 states of the bodies to the conventional form.
    /*!
     * Function to convert the USM6 elements state to the conventional form. For the USM6
     * propagator, this transforms USM6 elements w.r.t. the central bodies to the Cartesian states w.r.t. these
     * same central bodies: In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in USM6 elemements (i.e. form that is used in
     * numerical integration)
     * \param time Current time at which the state is valid
     * \param currentCartesianLocalSoluton State (internalSolution, which is Encke-formulation),
     *  converted to the 'conventional form' (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        // Convert state to Cartesian for each body
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            currentCartesianLocalSoluton.segment( i * 6, 6 ) =
                    orbital_element_conversions::convertUnifiedStateModelModifiedRodriguesParametersToCartesianElements(
                        internalSolution.block( i * 7, 0, 7, 1 ).template cast< double >( ), static_cast< double >(
                            centralBodyGravitationalParameters_.at( i )( ) ) ).template cast< StateScalarType >( );
        }

        currentCartesianLocalSoluton_ = currentCartesianLocalSoluton;
    }

    //! Function to get the acceleration models
    /*!
     * Function to get the acceleration models, including the central body accelerations that are removed for the
     * propagation scheme
     * \return List of acceleration models, including the central body accelerations that are removed in this propagation scheme.
     */
    basic_astrodynamics::AccelerationMap getFullAccelerationsMap( )
    {
        return originalAccelerationModelsPerBody_;
    }

    int getPropagatedStateSize( )
    {
        return 7 * this->bodiesToBeIntegratedNumerically_.size( );
    }

    void postProcessState(
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& unnormalizedState,
            const int startRow )
    {
        // Loop over each body
        for( unsigned int i = 0; i < this->bodiesToBeIntegratedNumerically_.size( ); i++ )
        {
            // Convert to/from shadow modifed Rodrigues parameters (SMRP) (transformation is the same either way)
            Eigen::Matrix< StateScalarType, 3, 1 > modifiedRodriguesParametersVector =
                    unnormalizedState.block( startRow + i * 7 + 3, 0, 3, 1 );
            StateScalarType modifiedRodriguesParametersMagnitude = modifiedRodriguesParametersVector.norm( );
            if ( modifiedRodriguesParametersMagnitude >= 1.0 )
            {
                // Invert flag
                unnormalizedState( startRow + i * 7 + 6 ) = std::fabs( unnormalizedState( startRow + i * 7 + 6 ) - 1.0 );

                // Convert to MRP/SMRP
                modifiedRodriguesParametersVector /= - modifiedRodriguesParametersMagnitude *
                        modifiedRodriguesParametersMagnitude;

                // Replace MRP with SMPR, or vice-versa
                unnormalizedState.segment( startRow + i * 7 + 3, 3 ) = modifiedRodriguesParametersVector;
            }
        }
    }

    virtual bool isStateToBePostProcessed( )
    {
        return true;
    }


private:

    //!  Gravitational parameters of central bodies used to convert Cartesian to Keplerian orbits, and vice versa
    std::vector< boost::function< double( ) > > centralBodyGravitationalParameters_;

    //! Central body accelerations for each propagated body, which has been removed from accelerationModelsPerBody_
    std::vector< boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
    centralAccelerations_;

    //! List of acceleration models, including the central body accelerations that are removed in this propagation scheme.
    basic_astrodynamics::AccelerationMap originalAccelerationModelsPerBody_;

    //! Current full Cartesian state of the propagated bodies, w.r.t. the central bodies
    /*!
     *  Current full Cartesian state of the propagated bodies, w.r.t. the central bodies. These variables are set when calling
     *  the convertToOutputSolution function.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentCartesianLocalSoluton_;

};


} // namespace propagators

} // namespace tudat

#endif // TUDAT_NUNIFIEDSTATEMODELMODIFIEDRODRIGUESPARAMETERSSTATEDERIVATIVE_H