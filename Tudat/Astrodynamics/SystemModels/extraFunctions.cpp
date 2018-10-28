#include "Tudat/Astrodynamics/SystemModels/extraFunctions.h"

#include "Tudat/Basics/utilities.h"

#include "Tudat/Astrodynamics/Aerodynamics/customConstantTemperatureAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"
#include "Tudat/Mathematics/BasicMathematics/numericalDerivative.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace system_models
{

//! Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
double areaBisectionFunction( const double currentTimeGuess, const double constantTimeStep, const Eigen::VectorXd& onboardTime,
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

//! Function to be used as input to the root-finder to determine the magnitude of the apoapsis maneuver.
double maneuverBisectionFunction( const double currentMagnitudeGuess, const Eigen::Vector6d& initialEstimatedCartesianState,
                                  const double targetPeriapsisRadius, const Eigen::Matrix3d& transformationFromLocalToInertialFrame,
                                  const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
                                  std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction,
                                  const bool apoapsisMaenuverEstimation )
{
    // Copy initial Cartesian state
    Eigen::Vector6d initialCartesianState = initialEstimatedCartesianState;

    // Create apoapsis maneuver vector in local frame and convert to inertial frame
    Eigen::Vector3d apoapsisManeuver = Eigen::Vector3d::Zero( );
    apoapsisManeuver[ 1 ] = currentMagnitudeGuess;
    apoapsisManeuver = transformationFromLocalToInertialFrame * apoapsisManeuver;

    // Add maneuver to initial state
    initialCartesianState.segment( 3, 3 ) += apoapsisManeuver;

    // Propagate state for two thirds of the orbit
    std::map< double, Eigen::VectorXd > propagatedState = statePropagationFunction( initialCartesianState ).second.first;

    // Retrieve periapsis altitude
    unsigned int i = 0;
    Eigen::VectorXd historyOfRadialDistances;
    historyOfRadialDistances.resize( propagatedState.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = propagatedState.begin( );
          mapIterator != propagatedState.end( ); mapIterator++, i++ )
    {
        historyOfRadialDistances[ i ] = mapIterator->second.segment( 0, 3 ).norm( );
    }
    double predictedApsisRadius;
    if ( apoapsisMaenuverEstimation )
    {
        predictedApsisRadius = historyOfRadialDistances.minCoeff( );
    }
    else
    {
        predictedApsisRadius = historyOfRadialDistances.maxCoeff( );
    }

    // Return predicted apo- or periapsis altitude w.r.t. expected value
    return predictedApsisRadius - targetPeriapsisRadius;
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
        const Eigen::VectorXd& estimatedAltitudesBelowAtmosphericInterface,
        const double referenceAltitude )
{
    using mathematical_constants::PI;

    // Pre-allocate variables
    double relativeAltitude;
    double twoPiRelativeAltitude;
    Eigen::VectorXd expectedDensity;
    Eigen::MatrixXd jacobianMatrix;
    expectedDensity.resizeLike( estimatedAltitudesBelowAtmosphericInterface );
    jacobianMatrix.resize( estimatedAltitudesBelowAtmosphericInterface.rows( ), 5 );

    // Loop over each acceleration to add values to matrix
    for ( unsigned int i = 0; i < estimatedAltitudesBelowAtmosphericInterface.rows( ); i++ )
    {
        // Find current expected measurement
        expectedDensity[ i ] = std::log( aerodynamics::threeTermAtmosphereModel(
                                             estimatedAltitudesBelowAtmosphericInterface[ i ], 0.0, 0.0, 0.0,
                                             referenceAltitude, std::exp( currentParameterEstimate[ 0 ] ),
                                         1.0 / currentParameterEstimate[ 1 ],
                std::vector< double >{ currentParameterEstimate[ 2 ], currentParameterEstimate[ 3 ], currentParameterEstimate[ 4 ] } ) );

        // Find current Jacobian matrix
        relativeAltitude = estimatedAltitudesBelowAtmosphericInterface[ i ] - referenceAltitude;
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

//! Altitude correction function for guidance system corridor estimator.
double correctionFactorForCorridorBoundaries( const double periapsisCorridorAltitude, const Eigen::Vector2d& linearLeastSquaresEstimate )
{
    return std::max( 1.0, 1.0 / ( linearLeastSquaresEstimate[ 0 ] + linearLeastSquaresEstimate[ 1 ] * periapsisCorridorAltitude ) );
}

} // namespace system_models

} // namespace tudat
