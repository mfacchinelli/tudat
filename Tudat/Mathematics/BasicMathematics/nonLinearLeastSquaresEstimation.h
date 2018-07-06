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

#ifndef TUDAT_NONLINEAR_LEAST_SQUARES_ESTIMATION_H
#define TUDAT_NONLINEAR_LEAST_SQUARES_ESTIMATION_H

#include <map>

#include <Eigen/Core>
#include <boost/function.hpp>

namespace tudat
{

namespace linear_algebra
{

//! Function to perform a non-linear least squares estimation.
Eigen::VectorXd nonLinearLeastSquaresFit(
        const boost::function< std::pair< Eigen::VectorXd, Eigen::MatrixXd >( const Eigen::VectorXd& ) >& observationAndDesignMatrixFunctions,
        const Eigen::VectorXd& initialEstimate, const Eigen::VectorXd& actualObservations,
        const double convergenceTolerance = 1.0e-8, const unsigned int maximumNumberOfInerations = 10 );

} // namespace linear_algebra

} // namespace tudat

#endif // TUDAT_NONLINEAR_LEAST_SQUARES_ESTIMATION_H
