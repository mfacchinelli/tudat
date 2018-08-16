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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Fornberg, B., "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids",
 *          Mathematics of Computation, October 1988.
 *
 */

#include "Tudat/Mathematics/BasicMathematics/numericalDerivative.h"

namespace tudat
{

namespace numerical_derivatives
{

//! Get coefficients of a certain order for central difference numerical derivatives.
std::map< int, double > getCentralDifferenceCoefficients( const CentralDifferenceOrders order )
{
    // Give value to the coefficients map, based on order
    std::map< int, double > coefficients;
    switch ( order )
    {
    case order2:
    {
        coefficients[ -1 ] = -1.0 / 2.0;
        coefficients[ 1 ] = 1.0 / 2.0;
        break;
    }
    case order4:
    {
        coefficients[ -2 ] = 1.0 / 12.0;
        coefficients[ -1 ] = -2.0 / 3.0;
        coefficients[ 1 ] = 2.0 / 3.0;
        coefficients[ 2 ] = -1.0 / 12.0;
        break;
    }
    case order6:
    {
        coefficients[ -3 ] = -1.0 / 60.0;
        coefficients[ -2 ] = 3.0 / 20.0;
        coefficients[ -1 ] = -3.0 / 4.0;
        coefficients[ 1 ] = 1.0 / 60.0;
        coefficients[ 2 ] = -3.0 / 20.0;
        coefficients[ 3 ] = 3.0 / 4.0;
        break;
    }
    case order8:
    {
        coefficients[ -4 ] = 1.0 / 280.0;
        coefficients[ -3 ] = -4.0 / 105.0;
        coefficients[ -2 ] = 1.0 / 5.0;
        coefficients[ -1 ] = -4.0 / 5.0;
        coefficients[ 1 ] = 4.0 / 5.0;
        coefficients[ 2 ] = -1.0 / 5.0;
        coefficients[ 3 ] = 4.0 / 105.0;
        coefficients[ 4 ] = -1.0 / 280.0;
        break;
    }
    }
    return coefficients;
}

Eigen::MatrixXd computeCentralDifference( const Eigen::VectorXd& input, const boost::function<
                                          Eigen::VectorXd( const Eigen::VectorXd& ) >& function,
                                          double minimumStep, double relativeStepSize,
                                          const CentralDifferenceOrders order )
{
    Eigen::MatrixXd result;

    for ( int derivative = 0; derivative < input.rows( ); derivative++ )
    {
        Eigen::VectorXd partial = computeCentralDifference( input, derivative, function,
                                                            minimumStep, relativeStepSize, order );
        if ( result.size( ) == 0 )
        {
            result = Eigen::MatrixXd( partial.rows( ), input.rows( ) );
        }
        result.col( derivative ) = partial;
    }
    return result;
}

//! Function to compute central difference with double as output and Eigen::VectorXd as input.
double computeCentralDifference(
        const boost::function< double( const Eigen::VectorXd& ) >& dependentVariableFunction,
        const Eigen::VectorXd& nominalIndependentVariable,
        const Eigen::VectorXd& independentVariableStepSize,
        const CentralDifferenceOrders order )
{
    // Retrieve coefficients
    const std::map< int, double > coefficients = getCentralDifferenceCoefficients( order );

    // Pre-allocate variables
    Eigen::VectorXd perturbedInput;
    double perturbedOutput;
    double numericalDerivative = 0.0;

    // Compute the numerical derivative.
    for ( std::map< int, double >::const_iterator coefficientIterator = coefficients.begin( );
          coefficientIterator != coefficients.end( ); coefficientIterator++ )
    {
        // Generate deviated input.
        perturbedInput = nominalIndependentVariable + coefficientIterator->first * independentVariableStepSize;
        perturbedOutput = dependentVariableFunction( perturbedInput );

        // Compute derivative
        numericalDerivative += perturbedOutput * ( coefficientIterator->second / independentVariableStepSize.norm( ) );
    }
    return numericalDerivative;
}

} // namespace numerical_derivatives

} // namespace tudat


