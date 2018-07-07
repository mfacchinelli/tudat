#include "Tudat/Astrodynamics/GuidanceNavigationControl/extraFunctions.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"
#include "Tudat/Astrodynamics/Aerodynamics/customConstantTemperatureAtmosphere.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
double areaBisectionFunction( const double currentTimeGuess, const double constantTimeStep,
                              const Eigen::VectorXd& onboardTime,
                              const std::vector< double >& estimatedAerodynamicAccelerationMagnitude )
{
    // Find nearest lower index to true anomaly guess
    int nearestLowerIndex = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch( onboardTime, currentTimeGuess );
    nearestLowerIndex = ( nearestLowerIndex == 0 ) ? 1 : nearestLowerIndex;

    // Compute trapezoidal quadrature to integrate the area until and after the current guess
    double lowerSliceQuadratureResult = numerical_quadrature::performExtendedSimpsonsQuadrature(
                constantTimeStep, utilities::sliceStlVector( estimatedAerodynamicAccelerationMagnitude, 0, nearestLowerIndex ) );
    double upperSliceQuadratureResult = numerical_quadrature::performExtendedSimpsonsQuadrature(
                constantTimeStep, utilities::sliceStlVector( estimatedAerodynamicAccelerationMagnitude, nearestLowerIndex ) );

    // Return difference in areas
    return upperSliceQuadratureResult - lowerSliceQuadratureResult;
}

//! Function to be used as input to the non-linear least squares process to determine the accelerometer errors.
std::pair< Eigen::VectorXd, Eigen::MatrixXd > accelerometerErrorEstimationFunction(
        const Eigen::Vector6d& currentErrorEstimate,
        const std::vector< Eigen::Vector3d >& vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface )
{
    // Set variables to zero
    Eigen::VectorXd expectedAcceleration = Eigen::VectorXd::Zero(
                3 * vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ) );
    Eigen::MatrixXd jacobianMatrix = Eigen::MatrixXd::Zero(
                3 * vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ), 6 );

    // Loop over each acceleration to add values to matrix
    for ( unsigned int i = 0; i < vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.size( ); i++ )
    {
        // Find current expected measurement
        expectedAcceleration.segment( 3 * i, 3 ) =
                ( Eigen::Matrix3d::Identity( ) - Eigen::Matrix3d( currentErrorEstimate.segment( 3, 3 ).asDiagonal( ) ) ) *
                ( vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i ) - currentErrorEstimate.segment( 0, 3 ) );

        // Find current Jacobian matrix
        jacobianMatrix.row( 3 * i ) << currentErrorEstimate[ 3 ] - 1.0, 0.0, 0.0, currentErrorEstimate[ 0 ] -
                vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i )[ 0 ], 0.0, 0.0;
        jacobianMatrix.row( 3 * i + 1 ) << 0.0, currentErrorEstimate[ 4 ] - 1.0, 0.0, 0.0, currentErrorEstimate[ 1 ] -
                vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i )[ 1 ], 0.0;
        jacobianMatrix.row( 3 * i + 2 ) << 0.0, 0.0, currentErrorEstimate[ 5 ] - 1.0, 0.0, 0.0, currentErrorEstimate[ 2 ] -
                vectorOfMeasuredAerodynamicAccelerationBelowAtmosphericInterface.at( i )[ 2 ];
    }

    // Return acceleration and design matrix as pair
    return std::make_pair( expectedAcceleration, jacobianMatrix );
}

//! Function to be used as input to the non-linear least squares process to determine the atmosphere parameters of the
//! three-term atmosphere model.
std::pair< Eigen::VectorXd, Eigen::MatrixXd > threeModelParametersEstimationFunction(
        const Eigen::Vector5d& currentParameterEstimate,
        const std::vector< double >& vectorOfEstimatedAltitudesBelowAtmosphericInterface,
        const double referenceAltitude )
{
    using mathematical_constants::PI;
    std::cout << "Parameters: " << currentParameterEstimate.transpose( ) << std::endl;

    // Set variables to zero
    Eigen::VectorXd expectedDensity = Eigen::VectorXd::Zero( vectorOfEstimatedAltitudesBelowAtmosphericInterface.size( ) );
    Eigen::MatrixXd jacobianMatrix = Eigen::MatrixXd::Zero( vectorOfEstimatedAltitudesBelowAtmosphericInterface.size( ), 5 );

    // Pre-allocate variables
    double relativeAltitude;
    double twoPiRelativeAltitude;

    // Loop over each acceleration to add values to matrix
    for ( unsigned int i = 0; i < vectorOfEstimatedAltitudesBelowAtmosphericInterface.size( ); i++ )
    {
        // Find current expected measurement
        expectedDensity[ i ] = std::log( aerodynamics::threeTermAtmosphereModel(
                    vectorOfEstimatedAltitudesBelowAtmosphericInterface.at( i ), 0.0, 0.0, 0.0, std::exp( currentParameterEstimate[ 0 ] ),
                referenceAltitude, 1.0 / currentParameterEstimate[ 1 ], { currentParameterEstimate[ 2 ],
                    currentParameterEstimate[ 3 ], currentParameterEstimate[ 4 ] } ) );

        // Find current Jacobian matrix
        relativeAltitude = vectorOfEstimatedAltitudesBelowAtmosphericInterface.at( i ) - referenceAltitude;
        twoPiRelativeAltitude = 2.0 * PI * relativeAltitude;
        jacobianMatrix.row( i ) << 1.0,
                currentParameterEstimate[ 2 ] * relativeAltitude -
                currentParameterEstimate[ 3 ] * std::sin( twoPiRelativeAltitude * currentParameterEstimate[ 1 ] ) * twoPiRelativeAltitude +
                currentParameterEstimate[ 4 ] * std::cos( twoPiRelativeAltitude * currentParameterEstimate[ 1 ] ) * twoPiRelativeAltitude,
                currentParameterEstimate[ 1 ] * relativeAltitude,
                std::cos( twoPiRelativeAltitude * currentParameterEstimate[ 1 ] ),
                std::sin( twoPiRelativeAltitude * currentParameterEstimate[ 1 ] );
    }

    // Return acceleration and design matrix as pair
    return std::make_pair( expectedDensity, jacobianMatrix );
}

} // namespace navigation

} // namespace tudat
