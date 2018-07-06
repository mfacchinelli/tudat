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

#include <cmath>
#include <iostream>

#include <Eigen/LU>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Mathematics/BasicMathematics/nonLinearLeastSquaresEstimation.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to perform a non-linear least squares estimation.
Eigen::VectorXd nonLinearLeastSquaresFit(
        const boost::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndDesignMatrixFunctions,
        const Eigen::VectorXd& initialEstimate, const Eigen::VectorXd& actualObservations,
        const double convergenceTolerance, const unsigned int maximumNumberOfInerations )
{
    // Set current estimate to initial value
    Eigen::VectorXd currentEstimate = initialEstimate;

    // Initialize variables
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > pairOfEstimatedObservationsAndDesignMatrix;
    Eigen::MatrixXd designMatrix;
    Eigen::VectorXd offsetInObservations;
    Eigen::VectorXd updateInEstimate;

    std::cout << "Starting iterative process." << std::endl;

    // Start iterative loop
    unsigned int iteration = 0;
    do
    {
        // Compute current system and jacobian functions
        pairOfEstimatedObservationsAndDesignMatrix = observationAndDesignMatrixFunctions( currentEstimate );
        designMatrix = pairOfEstimatedObservationsAndDesignMatrix.second;

        // Offset in observation
        offsetInObservations = actualObservations - pairOfEstimatedObservationsAndDesignMatrix.first;

        // Compute update in estimate
        updateInEstimate = ( designMatrix.transpose( ) * designMatrix ).inverse( ) * designMatrix.transpose( ) * offsetInObservations;
        std::cout << "Iteration: " << iteration << ". Update: " << updateInEstimate.transpose( ) << std::endl;

        // Correct estimate
        currentEstimate += updateInEstimate;
        iteration++;
    }
    while ( ( updateInEstimate.norm( ) > convergenceTolerance ) || ( iteration < maximumNumberOfInerations ) );

    // Give out new estimate in parameters
    std::cout << "Final: " << currentEstimate.transpose( ) << std::endl;
    return currentEstimate;
}

} // namespace linear_algebra

} // namespace tudat
