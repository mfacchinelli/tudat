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

#ifndef TUDAT_FILTER_H
#define TUDAT_FILTER_H

#include <map>
#include <limits>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

#include "Tudat/Basics/identityElements.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace filters
{

//! Filter class.
/*!
 *  Base class for the set up and use of filters.
 *  \tparam IndependentVariableType Type of independent variable. Default is double.
 *  \tparam DependentVariableType Type of dependent variable. Default is double.
 */
template< typename IndependentVariableType = double, typename DependentVariableType = double >
class FilterBase
{
public:

    //! Typedef of the state and measurement vectors.
    typedef Eigen::Matrix< DependentVariableType, Eigen::Dynamic, 1 > DependentVector;

    //! Typedef of the state and measurement matrices.
    typedef Eigen::Matrix< DependentVariableType, Eigen::Dynamic, Eigen::Dynamic > DependentMatrix;

    //! Typedef of the function describing the system and the measurements.
    typedef boost::function< DependentVector( const IndependentVariableType,
                                              const DependentVector& ) > Function;

    //! Typedefs for system and measurement matrix functions.
    typedef boost::function< DependentMatrix( const IndependentVariableType,
                                              const DependentVector& ) > MatrixFunction;

    //! Typedef of the integrator settings.
    typedef numerical_integrators::IntegratorSettings< IndependentVariableType > IntegratorSettings;

    //! Typedef of the integrator.
    typedef numerical_integrators::NumericalIntegrator< IndependentVariableType, DependentVector,
    DependentVector, IndependentVariableType > Integrator;

    //! Constructor.
    /*!
     *  Constructor.
     *  \param systemUncertainty Matrix defining the uncertainty in modeling of the system.
     *  \param measurementUncertainty Matrix defining the uncertainty in modeling of the measurements.
     *  \param initialTime Scalar representing the value of the initial time.
     *  \param initialStateVector Vector representing the initial (estimated) state of the system. It is used as first
     *      a-priori estimate of the state vector.
     *  \param initialCovarianceMatrix Matrix representing the initial (estimated) covariance of the system. It is used as first
     *      a-priori estimate of the covariance matrix.
     *  \param isStateToBeIntegrated Boolean defining whether the system function needs to be integrated.
     *  \param integrator Pointer to integrator to be used to propagate state.
     */
    FilterBase( const DependentMatrix& systemUncertainty,
                const DependentMatrix& measurementUncertainty,
                const IndependentVariableType initialTime,
                const DependentVector& initialStateVector,
                const DependentMatrix& initialCovarianceMatrix,
                const boost::shared_ptr< IntegratorSettings > integratorSettings ) :
        systemUncertainty_( systemUncertainty ), measurementUncertainty_( measurementUncertainty ), initialTime_( initialTime ),
        aPosterioriStateEstimate_( initialStateVector ), aPosterioriCovarianceEstimate_( initialCovarianceMatrix )
    {
        // Check that uncertainty matrices are square
        if ( systemUncertainty_.rows( ) != systemUncertainty_.cols( ) )
        {
            throw std::runtime_error( "Error in setting up filter. The system uncertainty matrix has to be square." );
        }
        if ( measurementUncertainty_.rows( ) != measurementUncertainty_.cols( ) )
        {
            throw std::runtime_error( "Error in setting up filter. The measurement uncertainty matrix has to be square." );
        }

        // Create noise distributions
        generateNoiseDistributions( );

        // Create system and measurement functions based on input parameters
        systemFunction_ = boost::bind( &FilterBase< IndependentVariableType,
                                       DependentVariableType >::createSystemFunction,
                                       this, _1, _2 );
        measurementFunction_ = boost::bind( &FilterBase< IndependentVariableType,
                                            DependentVariableType >::createMeasurementFunction,
                                            this, _1, _2 );

        // Create numerical integrator
        isStateToBeIntegrated_ = integratorSettings != nullptr;
        if ( isStateToBeIntegrated_ )
        {
            generateNumericalIntegrator( integratorSettings );
        }

        // Generate identity matrix
        identityMatrix_ = DependentMatrix::Identity( systemUncertainty_.rows( ), systemUncertainty_.cols( ) );

        // Add initial values to history
        historyOfStateEstimates_[ initialTime ] = aPosterioriStateEstimate_;
        historyOfCovarianceEstimates_[ initialTime ] = aPosterioriCovarianceEstimate_;
    }

    //! Destructor.
    virtual ~FilterBase( ){ }

    //! Function to update the filter with the data from the new time step.
    /*!
     *  Function to update the filter with the new step data.
     *  \param currentTime Scalar representing current time.
     *  \param currentMeasurementVector Vector representing current measurement.
     */
    virtual void updateFilter( const IndependentVariableType currentTime,
                               const DependentVector& currentMeasurementVector ) = 0;

    //! Function to update the a-posteriori estimates of state and covariance with external data.
    void modifyCurrentStateAndCovarianceEstimates( const DependentVector& newStateEstimate,
                                                   const DependentMatrix& newCovarianceEstimate )
    {
        // Update estimates with user-provided data
        aPosterioriStateEstimate_ = newStateEstimate;
        aPosterioriCovarianceEstimate_ = newCovarianceEstimate;
    }

    //! Function to produce system noise.
    /*!
     *  Function to produce system noise, based on a Gaussian distribution, with zero mean and standard
     *  deviation given by the diagonal elements of the input system uncertainty matrix.
     *  \return Vector representing system noise.
     */
    DependentVector produceSystemNoise( )
    {
        // Declare system noise vector
        DependentVector systemNoise = DependentVector::Zero( systemUncertainty_.rows( ) );

        // Loop over dimensions and add noise
        for ( int i = 0; i < systemUncertainty_.rows( ); i++ )
        {
            if ( systemNoiseDistribution_.at( i ) != nullptr )
            {
                systemNoise[ i ] = static_cast< DependentVariableType >(
                            systemNoiseDistribution_.at( i )->getRandomVariableValue( ) );
            }
        }

        // Give back noise
        systemNoiseHistory_.push_back( systemNoise );
        return systemNoise;
    }

    //! Function to produce measurement noise.
    /*!
     *  Function to produce measurement noise, based on a Gaussian distribution, with zero mean and standard
     *  deviation given by the diagonal elements of the input measurement uncertainty matrix.
     *  \return Vector representing measurement noise.
     */
    DependentVector produceMeasurementNoise( )
    {
        // Declare measurement noise vector
        DependentVector measurementNoise = DependentVector::Zero( measurementUncertainty_.rows( ) );

        // Loop over dimensions and add noise
        for ( int i = 0; i < measurementUncertainty_.rows( ); i++ )
        {
            if ( measurementNoiseDistribution_.at( i ) != nullptr )
            {
                measurementNoise[ i ] = static_cast< DependentVariableType >(
                            measurementNoiseDistribution_.at( i )->getRandomVariableValue( ) );
            }
        }

        // Give back noise
        measurementNoiseHistory_.push_back( measurementNoise );
        return measurementNoise;
    }

    //! Function to retrieve initial time.
    IndependentVariableType getInitialTime( ) { return initialTime_; }

    //! Function to retrieve step size for integration.
    IndependentVariableType getIntegrationStepSize( ) { return integrationStepSize_; }

    //! Function to retrieve current state estimate.
    /*!
     *  Function to retrieve current state estimate. The state estimate needs to first be computed by updating the
     *  filter with the updateFilter function.
     *  \return Current state estimate.
     */
    DependentVector getCurrentStateEstimate( ) { return aPosterioriStateEstimate_; }

    //! Function to retrieve current covariance estimate.
    /*!
     *  Function to retrieve current covariance estimate. The covariance estimate needs to first be computed by
     *  updating the filter with the updateFilter function.
     *  \return Current state estimate.
     */
    DependentMatrix getCurrentCovarianceEstimate( ) { return aPosterioriCovarianceEstimate_; }

    //! Function to retrieve the history of estimated states.
    /*!
     *  Function to retrieve the history of estimated states.
     *  \return History of estimated states for each time step.
     */
    std::map< IndependentVariableType, DependentVector > getEstimatedStateHistory( )
    {
        return historyOfStateEstimates_;
    }

    //! Function to retrieve the history of estimated covariance matrices.
    /*!
     *  Function to retrieve the history of estimated covariance matrices.
     *  \return History of estimated covariance matrices for each time step.
     */
    std::map< IndependentVariableType, DependentMatrix > getEstimatedCovarianceHistory( )
    {
        return historyOfCovarianceEstimates_;
    }

    //! Function to retrieve the history of system and measurement noise used by the updateFilter function.
    /*!
     *  Function to retrieve the history of system and measurement noise.
     *  \return History of system and measurement noise for each time step, output as a std::pair.
     */
    std::pair< std::vector< DependentVector >, std::vector< DependentVector > > getNoiseHistory( )
    {
        return std::make_pair( systemNoiseHistory_, measurementNoiseHistory_ );
    }

    //! Function to clear the history of stored variables.
    /*!
     *  Function to clear the history of stored variables. This function should be called if the history of state and covariance
     *  estimates over time needs to be deleted. This may be useful in case the filter is run for very long times.
     */
    virtual void clearFilterHistory( )
    {
        historyOfStateEstimates_.clear( );
        historyOfCovarianceEstimates_.clear( );
    }

protected:

    //! Function to create the function that defines the system model.
    /*!
     *  Function to create the function that defines the system model. The output of this function is then bound
     *  to the systemFunction_ variable, via the boost::bind command.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state.
     *  \return Vector representing the estimated state.
     */
    virtual DependentVector createSystemFunction( const IndependentVariableType currentTime,
                                                  const DependentVector& currentStateVector ) = 0;

    //! Function to create the function that defines the measurement model.
    /*!
     *  Function to create the function that defines the measurement model. The output of this function is then bound
     *  to the measurementFunction_ variable, via the boost::bind command.
     *  \param currentTime Scalar representing the current time.
     *  \param currentStateVector Vector representing the current state.
     *  \return Vector representing the estimated measurement.
     */
    virtual DependentVector createMeasurementFunction( const IndependentVariableType currentTime,
                                                       const DependentVector& currentStateVector ) = 0;

    //! Function to predict the state for the next time step.
    /*!
     *  Function to predict the state for the next time step, with the either the use of the integrator provided in
     *  the integratorSettings, or the systemFunction_ input by the user.
     *  \param currentTime Scalar representing the current time.
     *  \return Propagated state at the requested time.
     */
    virtual DependentVector predictState( const IndependentVariableType currentTime ) = 0;

    //! Function to correct the state for the next time step.
    /*!
     *  Function to predict the state for the next time step, by overwriting previous state, with the either the use of
     *  the integrator provided in the integratorSettings, or the systemFunction_ input by the user.
     *  \param currentTime Scalar representing the current time.
     *  \param aPrioriStateEstimate Vector denoting the a-priori state estimate.
     *  \param currentMeasurementVector Vector denoting the external measurement.
     *  \param measurementEstimate Vector denoting the measurement estimate.
     *  \param gainMatrix Gain matrix, such as Kalman gain (for Kalman filters).
     */
    void correctState( const IndependentVariableType currentTime,
                       const DependentVector& aPrioriStateEstimate, const DependentVector& currentMeasurementVector,
                       const DependentVector& measurementEstimate, const DependentMatrix& gainMatrix )
    {
        aPosterioriStateEstimate_ = aPrioriStateEstimate + gainMatrix * ( currentMeasurementVector - measurementEstimate );
        historyOfStateEstimates_[ currentTime ] = aPosterioriStateEstimate_;
    }

    //! Function to correct the covariance for the next time step.
    /*!
     *  Function to predict the state for the next time step, by overwriting previous state, with the either the use of
     *  the integrator provided in the integratorSettings, or the systemFunction_ input by the user.
     *  \param currentTime Scalar representing the current time.
     */
    virtual void correctCovariance( const IndependentVariableType currentTime, const DependentMatrix& aPrioriCovarianceEstimate,
                                    const DependentMatrix& currentMeasurementMatrix, const DependentMatrix& kalmanGain ) = 0;

    //! System function.
    /*!
     *  System function that will be used to retrieve the a-priori estimated state for the next step.
     */
    Function systemFunction_;

    //! Measurement function.
    /*!
     *  Measurement function that will be used to retrieve the estimated measurement for the next step,
     *  based on the current state.
     */
    Function measurementFunction_;

    //! Matrix representing the uncertainty in system modeling.
    DependentMatrix systemUncertainty_;

    //! Matrix representing the uncertainty in measurement modeling.
    DependentMatrix measurementUncertainty_;

    //! Scalar representing the initial time.
    IndependentVariableType initialTime_;

    //! Vector representing the a-posteriori estimated state.
    /*!
     *  Vector representing the a-posteriori estimated state, i.e., the state after the prediction and
     *  update steps of the Kalman filter.
     */
    DependentVector aPosterioriStateEstimate_;

    //! Matrix representing the a-posteriori estimated covariance.
    /*!
     *  Matrix representing the a-posteriori estimated covariance, i.e., the covariance after the prediction and
     *  update steps of the Kalman filter.
     */
    DependentMatrix aPosterioriCovarianceEstimate_;

    //! Boolean specifying whether the state needs to be integrated.
    bool isStateToBeIntegrated_;

    //! Pointer to the integrator.
    /*!
     *  Pointer to the integrator, which is used to propagate the state to the new time step.
     */
    boost::shared_ptr< Integrator > integrator_;

    //! Scalar representing step size for integration.
    /*!
     *  Scalar representing step size for integration. If integrator_ points to a constant step size integrator, then
     *  this will be the constant step size, otherwise it will be the initial step size.
     */
    IndependentVariableType integrationStepSize_;

    //! Indentity matrix.
    /*!
     *  Indentity matrix with the correct dimensions for the specific application.
     */
    DependentMatrix identityMatrix_;

    //! Map of estimated states vectors history.
    std::map< IndependentVariableType, DependentVector > historyOfStateEstimates_;

    //! Map of estimated covariance matrices history.
    std::map< IndependentVariableType, DependentMatrix > historyOfCovarianceEstimates_;

private:

    //! Function to generate the noise distributions for both system and measurement modeling.
    /*!
     *  Function to generate the noise distributions for both system and measurement modeling, which uses
     *  a Gaussian distribution, with zero mean and standard deviation given by the diagonal elements of the
     *  input system and measurement uncertainty matrices.
     */
    void generateNoiseDistributions( )
    {
        using namespace tudat::statistics;

        // Create system noise
        for ( unsigned int i = 0; i < systemUncertainty_.rows( ); i++ )
        {
            if ( static_cast< double >( systemUncertainty_( i, i ) ) != 0.0 )
            {
                systemNoiseDistribution_.push_back(
                            createBoostContinuousRandomVariableGenerator(
                                normal_boost_distribution, { 0.0, static_cast< double >(
                                                             std::sqrt( systemUncertainty_( i, i ) ) ) },
                                static_cast< double >( i ) ) );
            }
            else
            {
                systemNoiseDistribution_.push_back( nullptr );
            }
        }

        // Create measurement noise
        for ( unsigned int i = 0; i < measurementUncertainty_.rows( ); i++ )
        {
            if ( static_cast< double >( measurementUncertainty_( i, i ) ) != 0.0 )
            {
                measurementNoiseDistribution_.push_back(
                            createBoostContinuousRandomVariableGenerator(
                                normal_boost_distribution, { 0.0, static_cast< double >(
                                                             std::sqrt( measurementUncertainty_( i, i ) ) ) },
                                static_cast< double >( systemUncertainty_.rows( ) + i ) ) );
            }
            else
            {
                measurementNoiseDistribution_.push_back( nullptr );
            }
        }
    }

    //! Function to generate the numerical integrator to be used for propagation of the state.
    /*!
     *  Function to generate the numerical integrator to be used for propagation of the state, based on the integrator
     *  settings and the systemFunction_ input by the user. The systemFunction_ therefore acts as the differential equation
     *  for the system. Currently, only Euler integration is supported.
     *  \param integratorSettings Pointer to integration settings.
     */
    void generateNumericalIntegrator( const boost::shared_ptr< IntegratorSettings > integratorSettings )
    {
        // Get time step information
        integrationStepSize_ = integratorSettings->initialTimeStep_;

        // Generate integrator
        switch ( integratorSettings->integratorType_ )
        {
        case numerical_integrators::euler:
        case numerical_integrators::rungeKutta4:
        {
            integrator_ = numerical_integrators::createIntegrator< IndependentVariableType, DependentVector >(
                        systemFunction_, aPosterioriStateEstimate_, integratorSettings );
            break;
        }
        default:
            throw std::runtime_error( "Error in setting up filter. Only constant time step integrators (i.e., Euler and "
                                      "Runge-Kutta 4) are supported." );
        }
    }

    //! Vector where the system noise generators are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > systemNoiseDistribution_;

    //! Vector where the measurement noise generators are stored.
    std::vector< boost::shared_ptr< statistics::RandomVariableGenerator< double > > > measurementNoiseDistribution_;

    //! Vector of system noise hisotries.
    std::vector< DependentVector > systemNoiseHistory_;

    //! Vector of measurement noise hisotries.
    std::vector< DependentVector > measurementNoiseHistory_;

};

//! Typedef for a filter with double data type.
typedef FilterBase< > FilterDouble;

//! Typedef for a shared-pointer to a filter with double data type.
typedef boost::shared_ptr< FilterDouble > FilterDoublePointer;

} // namespace filters

} // namespace tudat

#endif // TUDAT_FILTER_H
