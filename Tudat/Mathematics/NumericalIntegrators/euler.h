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
 *
 */

#ifndef TUDAT_EULER_INTEGRATOR_H
#define TUDAT_EULER_INTEGRATOR_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalIntegrators/reinitializableNumericalIntegrator.h"

namespace tudat
{
namespace numerical_integrators
{

//! Class that implements the Euler integrator.
/*!
 * Class that implements the Euler, fixed order, fixed step size integrator.
 * \tparam StateType The type of the state. This type should support addition with
 *          StateDerivativeType
 * \tparam StateDerivativeType The type of the state derivative. This type should support
 *          multiplication with IndependentVariableType and doubles.
 * \tparam IndependentVariableType The type of the independent variable.
 * \sa NumericalIntegrator.
 */
template< typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd, typename TimeStepType = IndependentVariableType >
class EulerIntegrator :
        public numerical_integrators::ReinitializableNumericalIntegrator<
        IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
public:

    //! Typedef of the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef numerical_integrators::ReinitializableNumericalIntegrator<
    IndependentVariableType, StateType,
    StateDerivativeType, TimeStepType > ReinitializableNumericalIntegratorBase;

    //! Typedef to the state derivative function.
    /*!
     * Typedef to the state derivative function inherited from the base class.
     * \sa NumericalIntegrator::StateDerivativeFunction.
     */
    typedef typename ReinitializableNumericalIntegratorBase::NumericalIntegratorBase::
    StateDerivativeFunction StateDerivativeFunction;

    //! Default constructor.
    /*!
     * Default constructor, taking a state derivative function as argument.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    EulerIntegrator( const StateDerivativeFunction& stateDerivativeFunction,
                     const IndependentVariableType intervalStart,
                     const StateType& initialState )
        : ReinitializableNumericalIntegratorBase( stateDerivativeFunction ),
          currentIndependentVariable_( intervalStart ),
          currentState_( initialState ),
          lastIndependentVariable_( intervalStart )
    { }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step.
     * \return Step size to be used for the next step.
     */
    virtual TimeStepType getNextStepSize( ) const { return stepSize_; }

    //! Get current state.
    /*!
     * Returns the current state of the Euler integrator.
     * \return Current integrated state.
     */
    virtual StateType getCurrentState( ) const { return currentState_; }

    //! Returns the current independent variable.
    /*!
     * Returns the current value of the independent variable of the integrator.
     * \return Current independent variable.
     */
    virtual IndependentVariableType getCurrentIndependentVariable( ) const
    {
        return currentIndependentVariable_;
    }

    //! Perform a single Euler integration step.
    /*!
     * Performs a single Euler integration step using a step of size specified by stepSize. The
     * initial state for the step is internally set to the final state of the previous step. In
     * case this is the first step, the initial state is set to the initial state provided by the
     * user.
     * \param stepSize The size of the step to take.
     * \return The state at the end of the interval.
     */
    virtual StateType performIntegrationStep( const TimeStepType stepSize )
    {
        lastIndependentVariable_ = currentIndependentVariable_;
        lastState_ = currentState_;

        currentState_ += stepSize * this->stateDerivativeFunction_( currentIndependentVariable_, currentState_ );

        stepSize_ = stepSize;
        currentIndependentVariable_ += stepSize_;

        // Return the integration result.
        return currentState_;
    }

    //! Rollback the internal state to the last state.
    /*!
     * Performs rollback of the internal state to the last state. This function can only be called
     * once after calling integrateTo() or performIntegrationStep() unless specified otherwise by
     * implementations, and can not be called before any of these functions have been called. Will
     * return true if the rollback was successful, and false otherwise.
     * \return True if the rollback was successful.
     */
    virtual bool rollbackToPreviousState( )
    {
        if ( currentIndependentVariable_ == lastIndependentVariable_ )
        {
            return false;
        }

        currentIndependentVariable_ = lastIndependentVariable_;
        currentState_ = lastState_;
        return true;
    }

    //! Get previous independent variable.
    /*!
     * Returns the previoius value of the independent variable of the integrator.
     * \return Previous independent variable.
     */
    IndependentVariableType getPreviousIndependentVariable( )
    {
        return this->lastIndependentVariable_;
    }

    //! Get previous state value.
    /*!
     * Returns the previous value of the state.
     * \return Previous state
     */
    StateType getPreviousState( )
    {
        return this->lastState_;
    }

    //! Modify the state at the current value of the independent variable.
    /*!
     * Modify the state at the current value of the independent variable.
     * \param newState The new state to set the current state to.
     */
    void modifyCurrentState( const StateType& newState, const IndependentVariableType newTime = 0 )
    {
        this->currentState_ = newState;
        if ( newTime == 0 )
        {
            this->lastIndependentVariable_ = currentIndependentVariable_;
        }
        else
        {
            this->lastIndependentVariable_ = newTime;
        }
    }

protected:

    //! Last used step size.
    /*!
     * Last used step size, passed to either integrateTo() or performIntegrationStep().
     */
    TimeStepType stepSize_;

    //! Current independent variable.
    /*!
     * Current independent variable as computed by performIntegrationStep().
     */
    IndependentVariableType currentIndependentVariable_;

    //! Current state.
    /*!
     * Current state as computed by performIntegrationStep( ).
     */
    StateType currentState_;

    //! Last independent variable.
    /*!
     * Last independent variable value as computed by performIntegrationStep( ).
     */
    IndependentVariableType lastIndependentVariable_;

    //! Last state.
    /*!
     * Last state as computed by performIntegrationStep( ).
     */
    StateType lastState_;
};

//! Typedef of Euler integrator (state/state derivative = VectorXd, independent variable = double).
/*!
 * Typedef of an Euler integrator with VectorXds as state and state derivative and double as
 * independent variable.
 */
typedef EulerIntegrator< > EulerIntegratorXd;

//! Typedef of a scalar Euler integrator.
/*!
 * Typedef of a Euler integrator with doubles as state and state derivative and independent variable.
 */
typedef EulerIntegrator< double, double, double > EulerIntegratord;

//! Typedef for a shared-pointer to default Euler integrator.
/*!
 * Typedef for a shared-pointer to an Euler integrator with VectorXds as state and state derivative and double
 * as independent variable.
 */
typedef boost::shared_ptr< EulerIntegratorXd > EulerIntegratorXdPointer;

//! Typedef of pointer to a scalar Euler integrator.
/*!
 * Typedef of pointer to an Euler integrator with doubles as state and state derivative and
 * independent variable.
 */
typedef boost::shared_ptr< EulerIntegratord > EulerIntegratordPointer;

} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_EULER_INTEGRATOR_H
