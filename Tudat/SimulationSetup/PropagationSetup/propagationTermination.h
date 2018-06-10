/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H
#define TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H

#include <boost/shared_ptr.hpp>

#include "Tudat/SimulationSetup/PropagationSetup/propagationOutput.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"

namespace tudat
{

namespace propagators
{

//! Possible events that can trigger the termination of a propagation
enum PropagationTerminationReason
{
    propagation_never_run,
    unknown_propagation_termination_reason,
    termination_condition_reached,
    runtime_error_caught_in_propagation,
    nan_or_inf_detected_in_state
};

//! Base class for checking whether the numerical propagation is to be stopped at current time step or not
/*!
 *  Base class for checking whether the numerical propagation is to be stopped at current time step or not. Derived
 *  classes implement the various types of conditions (and associated threshold values) under which the propagation is to
 *  be stopped.
 */
class PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param terminationType Type of termination condition
     * \param terminateExactlyOnFinalCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    PropagationTerminationCondition(
            const PropagationTerminationTypes terminationType,
            const bool terminateExactlyOnFinalCondition = false ):
        terminationType_( terminationType ), terminateExactlyOnFinalCondition_( terminateExactlyOnFinalCondition ){ }

    //! Destructor
    virtual ~PropagationTerminationCondition( ){ }

    //! (Pure virtual) function to check whether the propagation should be stopped
    /*!
     * (Pure virtual) function to check whether the propagation should be stopped. Note that the accelerations and
     * environment must be updated (done automatically during numerical propagation) to check the stopping condition.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    virtual bool checkStopCondition( const double time, const double cpuTime ) = 0;

    //! Function to retrieve type of termination condition
    /*!
     *  Function to retrieve type of termination condition
     *  \return Type of termination condition
     */
    virtual PropagationTerminationTypes getTerminationType( )
    {
        return terminationType_;
    }

    //! Function to retrieve boolean to denote whether the propagation is to terminate exactly on the final condition
    /*!
     *  Function to retrieve boolean to denote whether the propagation is to terminate exactly on the final condition
     *  \return Boolean to denote whether the propagation is to terminate exactly on the final condition
     */
    bool getTerminateExactlyOnFinalCondition( )
    {
        return terminateExactlyOnFinalCondition_;
    }

protected:

    //! Type of termination condition
    PropagationTerminationTypes terminationType_;

    //! Boolean to denote whether the propagation is to terminate exactly on the final condition, or whether it is to terminate
    //! on the first step where it is violated.
    bool terminateExactlyOnFinalCondition_;

};

//! Class for stopping the propagation after a fixed amount of time (i.e. for certain independent variable value)
class FixedTimePropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stopTime Time at which the propagation is to stop.
     * \param propagationDirectionIsPositive Boolean denoting whether propagation is forward (if true) or backwards
     * (if false) in time.
     * \param terminateExactlyOnFinalCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    FixedTimePropagationTerminationCondition(
            const double stopTime,
            const bool propagationDirectionIsPositive,
            const bool terminateExactlyOnFinalCondition = false ):
        PropagationTerminationCondition( time_stopping_condition, terminateExactlyOnFinalCondition ),
        stopTime_( stopTime ),
        propagationDirectionIsPositive_( propagationDirectionIsPositive ){ }


    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the stopTime_ has been reached or not.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

    //! Function to retrieve time at which the propagation is to stop.
    /*!
     *  Function to retrieve time at which the propagation is to stop.
     *  \return Type of termination condition
     */
    double getStopTime( )
    {
        return stopTime_;
    }

private:

    //! Time at which the propagation is to stop.
    double stopTime_;

private:

    //!  Boolean denoting whether propagation is forward (if true) or backwards (if false) in time.
    bool propagationDirectionIsPositive_;
};

//! Class for stopping the propagation after a fixed amount of CPU time
class FixedCPUTimePropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param cpuStopTime CPU time at which the propagation is to stop.
     */
    FixedCPUTimePropagationTerminationCondition( const double cpuStopTime ) :
        PropagationTerminationCondition( cpu_time_stopping_condition, false ),
        cpuStopTime_( cpuStopTime ) { }


    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the stopTime_ has been reached or not.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

private:

    //! Time at which the propagation is to stop.
    double cpuStopTime_;

};

//! Class for stopping the propagation when a dependent variable reaches a given value (either upper or lower bound)
class SingleVariableLimitPropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dependentVariableSettings Settings for dependent variable that is to be checked
     * \param variableRetrievalFuntion Function returning the dependent variable.
     * \param limitingValue Value at which the propagation is to be stopped
     * \param useAsLowerBound Boolean denoting whether the propagation should stop if the dependent variable goes below
     * (if true) or above (if false) limitingValue
     * \param terminateExactlyOnFinalCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     * \param terminationRootFinderSettings Settings to create root finder used to converge on exact final condition.
     */
    SingleVariableLimitPropagationTerminationCondition(
            const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const boost::function< double( ) > variableRetrievalFuntion,
            const double limitingValue,
            const bool useAsLowerBound,
            const bool terminateExactlyOnFinalCondition = false,
            const boost::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings = NULL ):
        PropagationTerminationCondition(
            dependent_variable_stopping_condition, terminateExactlyOnFinalCondition ),
        dependentVariableSettings_( dependentVariableSettings ), variableRetrievalFuntion_( variableRetrievalFuntion ),
        limitingValue_( limitingValue ), useAsLowerBound_( useAsLowerBound ),
        terminationRootFinderSettings_( terminationRootFinderSettings )
    {
        if( ( terminateExactlyOnFinalCondition == false ) && ( terminationRootFinderSettings != NULL ) )
        {
            std::cerr<<"Warning, root finder provided to SingleVariableLimitPropagationTerminationCondition, but termination on final conditions set to false"<<std::endl;
        }
        if( ( terminateExactlyOnFinalCondition == true ) && doesRootFinderRequireDerivatives( terminationRootFinderSettings ) )
        {
            throw std::runtime_error( "Error when setting exact dependent variable termination, requested root finder requires derivatives; not available in state derivative model" );
        }
    }

    //! Destructor.
    ~SingleVariableLimitPropagationTerminationCondition( ){ }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. whether the given dependent variable has been
     * reached or not.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

    //! Function to return current difference between termination variable, and the value at which the propagation must terminate
    /*!
     * Function to return current difference between termination variable, and the value at which the propagation must terminate
     * \return Current difference between termination variable, and the value at which the propagation must terminate
     */
    double getStopConditionError( )
    {
        return variableRetrievalFuntion_( ) - limitingValue_;
    }

    //! Function to retrieve settings to create root finder used to converge on exact final condition.
    /*!
     *  Function to retrieve settings to create root finder used to converge on exact final condition.
     *  \return Settings to create root finder used to converge on exact final condition.
     */
    boost::shared_ptr< root_finders::RootFinderSettings > getTerminationRootFinderSettings( )
    {
        return terminationRootFinderSettings_;
    }

private:

    //! Settings for dependent variable that is to be checked
    boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings_;

    //! Function returning the dependent variable.
    boost::function< double( ) > variableRetrievalFuntion_;

    //! Value at which the propagation is to be stopped
    double limitingValue_;

    //! Boolean denoting whether the propagation should stop if the dependent variable goes below
    //! (if true) or above (if false) limitingValue
    bool useAsLowerBound_;

    //! Settings to create root finder used to converge on exact final condition.
    boost::shared_ptr< root_finders::RootFinderSettings > terminationRootFinderSettings_;
};

//! Class for stopping the propagation when one or all of a given set of stopping conditions is reached.
class HybridPropagationTerminationCondition: public PropagationTerminationCondition
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param propagationTerminationCondition List of termination conditions that are checked when calling
     * checkStopCondition is called.
     * \param fulFillSingleCondition Boolean denoting whether a single (if true) or all (if false) of the entries in the
     * propagationTerminationCondition_ should return true from the checkStopCondition function to stop the propagation
     * \param terminateExactlyOnFinalCondition Boolean to denote whether the propagation is to terminate exactly on the final
     * condition, or whether it is to terminate on the first step where it is violated.
     */
    HybridPropagationTerminationCondition(
            const std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition,
            const bool fulFillSingleCondition = 0,
            const bool terminateExactlyOnFinalCondition = 0 ):
        PropagationTerminationCondition( hybrid_stopping_condition, terminateExactlyOnFinalCondition ),
        propagationTerminationCondition_( propagationTerminationCondition ),
        fulFillSingleCondition_( fulFillSingleCondition )
    {
        isConditionMetWhenStopping_.resize( propagationTerminationCondition.size( ) );
    }

    //! Function to check whether the propagation is to be be stopped
    /*!
     * Function to check whether the propagation is to be be stopped, i.e. one or all (depending on value of
     * fulFillSingleCondition_) of the stopping conditions are fulfilled.
     * \param time Current time in propagation
     * \param cpuTime Current CPU time in propagation
     * \return True if propagation is to be stopped, false otherwise.
     */
    bool checkStopCondition( const double time, const double cpuTime );

    //! Function to retrieve list of termination conditions that are checked when calling checkStopCondition is called.
    /*!
     *  Function to retrieve list of termination conditions that are checked when calling checkStopCondition is called.
     *  \return List of termination conditions that are checked when calling checkStopCondition is called.
     */
    std::vector< boost::shared_ptr< PropagationTerminationCondition > > getPropagationTerminationConditions( )
    {
        return propagationTerminationCondition_;
    }

    //! Function to retrieve whether all or a single termination condition should be met
    /*!
     *  Function to retrieve whether all or a single termination condition should be met
     *  \return Boolean denoting whether a single (if true) or all (if false) of the entries in the
     *  propagationTerminationCondition_ should return true from the checkStopCondition function to stop the propagation.
     */
    bool getFulFillSingleCondition( )
    {
        return fulFillSingleCondition_;
    }

    std::vector< bool > getIsConditionMetWhenStopping( )
    {
        return isConditionMetWhenStopping_;
    }

private:

    //! List of termination conditions that are checked when calling checkStopCondition is called.
    std::vector< boost::shared_ptr< PropagationTerminationCondition > > propagationTerminationCondition_;

    //!  Boolean denoting whether a single (if true) or all (if false) of the entries in the propagationTerminationCondition_
    //!  should return true from the checkStopCondition function to stop the propagation.
    bool fulFillSingleCondition_;

    std::vector< bool > isConditionMetWhenStopping_;

};

//! Function to create propagation termination conditions from associated settings
/*!
 * Function to create propagation termination conditions from associated settings
 * \param terminationSettings Settings for propagation termination conditions
 * \param bodyMap List of body objects that contains all environment models
 * \param initialTimeStep Time step at first call of numerical integration.
 * \return Object used to check whether propagation is to be stopped or not.
 */
boost::shared_ptr< PropagationTerminationCondition > createPropagationTerminationConditions(
        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const double initialTimeStep );

//! Class for storing details on the propagation termination
class PropagationTerminationDetails
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param propagationTerminationReason Reason for termination
     * \param terminationOnExactCondition True if exact termination condition is used, false if not, -1 if neither is relevant
     */
    PropagationTerminationDetails( const PropagationTerminationReason propagationTerminationReason,
                                   const bool terminationOnExactCondition = -1 ):
        propagationTerminationReason_( propagationTerminationReason ),
        terminationOnExactCondition_( terminationOnExactCondition ){ }

    //! Function to retrieve reason for termination
    /*!
     * Function to retrieve reason for termination
     * \return Reason for termination
     */
    PropagationTerminationReason getPropagationTerminationReason( )
    {
        return propagationTerminationReason_;
    }

    //! Function to retrieve boolean to denote whether exact termination conditions are used.
    /*!
     * Function to retrieve boolean to denote whether exact termination conditions are used.
     * \return Boolean to denote whether exact termination conditions are used.
     */
    bool getTerminationOnExactCondition( )
    {
        return terminationOnExactCondition_;
    }

protected:

    //! Reason for termination
    PropagationTerminationReason propagationTerminationReason_;

    //! Boolean to denote whether exact termination conditions are used.
    /*!
     *  Boolean to denote whether exact termination conditions are used. True if exact termination condition is used,
     *  false if not, -1 if neither is relevant.
     */
    bool terminationOnExactCondition_;

};

//! Class for storing details on the propagation termination when using hybrid termination conditions
class PropagationTerminationDetailsFromHybridCondition: public PropagationTerminationDetails
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param terminationOnExactCondition True if exact termination condition is used, false if not.
     * \param terminationCondition Hybrid termination conditions that were used
     */
    PropagationTerminationDetailsFromHybridCondition(
            const bool terminationOnExactCondition,
            const boost::shared_ptr< HybridPropagationTerminationCondition > terminationCondition ):
        PropagationTerminationDetails( termination_condition_reached, terminationOnExactCondition ),
        isConditionMetWhenStopping_( terminationCondition->getIsConditionMetWhenStopping( ) ){ }

    //! Function to retrieve list of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
    /*!
     * Function to retrieve list of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
     * \return list of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
     */
    std::vector< bool > getWasConditionMetWhenStopping( )
    {
        if( terminationOnExactCondition_ )
        {
            std::cerr<<"Warning when retrieving list of conditions that were met in hybrid propagation termination details. Propagation was terminated on exact conditions using root finder, list of conditions may not be reliable"<<std::endl;
        }
        return isConditionMetWhenStopping_;
    }

private:

    //! List of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
    std::vector< bool > isConditionMetWhenStopping_;

};


} // namespace propagators

} // namespace tudat


#endif // TUDAT_PROPAGATIONTERMINATIONCONDITIONS_H
