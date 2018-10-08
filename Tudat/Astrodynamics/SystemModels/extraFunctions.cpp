#include "Tudat/Astrodynamics/SystemModels/extraFunctions.h"

#include "Tudat/Basics/utilities.h"

#include "Tudat/Astrodynamics/Aerodynamics/customConstantTemperatureAtmosphere.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"
#include "Tudat/Mathematics/BasicMathematics/numericalDerivative.h"

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

//! Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
double lowerAltitudeBisectionFunctionBasedOnHeatingConditions(
        const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
        const double planetaryRadius, const double planetaryGravitationalParameter,
        const double maximumHeatRate, const double maximumHeatLoad,
        const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
        std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
{
    // Copy initial Keplerian state
    Eigen::Vector6d initialKeplerianState = initialEstimatedKeplerianState;

    // Modify initial state to match new estimated periapsis altitude
    double estimatedApoapsisRadius = basic_astrodynamics::computeKeplerRadialDistance(
                initialKeplerianState[ 0 ], initialKeplerianState[ 1 ], initialKeplerianState[ 5 ] );
    double semiMajorAxis = 0.5 * ( estimatedApoapsisRadius + currentAltitudeGuess + planetaryRadius );
    double eccentricity = ( estimatedApoapsisRadius - currentAltitudeGuess - planetaryRadius ) /
            ( estimatedApoapsisRadius + currentAltitudeGuess + planetaryRadius );
    initialKeplerianState[ 0 ] = semiMajorAxis;
    initialKeplerianState[ 1 ] = eccentricity;

    // Propagate orbit to new condition and retrieve heating conditions
    std::map< double, Eigen::VectorXd > heatingConditions =
            statePropagationFunction( orbital_element_conversions::convertKeplerianToCartesianElements(
                                          initialKeplerianState, planetaryGravitationalParameter ) ).second.second;

    // Separate time and heat rate and find maximum heat rate
    unsigned int i = 0;
    Eigen::VectorXd historyOfHeatRates;
    historyOfHeatRates.resize( heatingConditions.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = heatingConditions.begin( );
          mapIterator != heatingConditions.end( ); mapIterator++, i++ )
    {
        historyOfHeatRates[ i ] = mapIterator->second[ 1 ];
    }
    double heatRate = historyOfHeatRates.maxCoeff( );

    // Compute heat load by integrating heat flux
    double heatLoad = numerical_quadrature::performExtendedSimpsonsQuadrature(
                10.0, utilities::convertEigenVectorToStlVector( historyOfHeatRates ) );

    // Compute offsets w.r.t. maximum allowed heat rate and heat load
    double offsetInHeatRate = heatRate - maximumHeatRate;
    double offsetInHeatLoad = heatLoad - maximumHeatLoad;

    // Return value to to indicate closeness to limiting value
    return ( offsetInHeatRate > offsetInHeatLoad ) ? offsetInHeatRate : offsetInHeatLoad;
}

//! Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
double lowerAltitudeBisectionFunctionBasedOnLifetimeCondition(
        const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
        const double planetaryRadius, const double planetaryGravitationalParameter, const double minimumLifetime,
        const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
        std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
{
    // Copy initial Keplerian state
    Eigen::Vector6d initialKeplerianState = initialEstimatedKeplerianState;

    // Modify initial state to match new estimated periapsis altitude
    double estimatedPeriapsisRadius = currentAltitudeGuess + planetaryRadius;
    double estimatedApoapsisRadius = basic_astrodynamics::computeKeplerRadialDistance(
                initialKeplerianState[ 0 ], initialKeplerianState[ 1 ], initialKeplerianState[ 5 ] );
    double semiMajorAxis = 0.5 * ( estimatedApoapsisRadius + estimatedPeriapsisRadius );
    double eccentricity = ( estimatedApoapsisRadius - estimatedPeriapsisRadius ) / ( estimatedApoapsisRadius + estimatedPeriapsisRadius );
    initialKeplerianState[ 0 ] = semiMajorAxis;
    initialKeplerianState[ 1 ] = eccentricity;

    // Propagate orbit to new condition and retrieve heating conditions
    std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
            std::map< double, Eigen::VectorXd > > > propagationResult = statePropagationFunction(
                orbital_element_conversions::convertKeplerianToCartesianElements(
                    initialKeplerianState, planetaryGravitationalParameter ) );

    // Give propagation a score based on lifetime
    double actualLifetime =
            ( propagationResult.second.first.rbegin( )->first - propagationResult.second.first.begin( )->first ) /
            physical_constants::JULIAN_DAY;
    std::cout << "Altitude: " << currentAltitudeGuess / 1.0e3 << "km. Lifetime: " << actualLifetime << " day" << std::endl
              << "Score: " << actualLifetime - minimumLifetime << std::endl << std::endl;

    // Return value to to indicate closeness to limiting value
    return actualLifetime - minimumLifetime;
}

//! Function to be used as input to the root-finder to determine the upper altitude bound for the periapsis corridor.
double upperAltitudeBisectionFunction( const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
                                       const double planetaryRadius, const double planetaryGravitationalParameter,
                                       const double minimumDynamicPressure,
                                       const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
                                       std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
{
    // Copy initial Keplerian state
    Eigen::Vector6d initialKeplerianState = initialEstimatedKeplerianState;

    // Modify initial state to match new estimated periapsis altitude
    double estimatedApoapsisRadius = basic_astrodynamics::computeKeplerRadialDistance(
                initialKeplerianState[ 0 ], initialKeplerianState[ 1 ], initialKeplerianState[ 5 ] );
    double semiMajorAxis = 0.5 * ( estimatedApoapsisRadius + currentAltitudeGuess + planetaryRadius );
    double eccentricity = ( estimatedApoapsisRadius - currentAltitudeGuess - planetaryRadius ) /
            ( estimatedApoapsisRadius + currentAltitudeGuess + planetaryRadius );
    initialKeplerianState[ 0 ] = semiMajorAxis;
    initialKeplerianState[ 1 ] = eccentricity;

    // Propagate orbit to new condition and retrieve heating conditions
    std::map< double, Eigen::VectorXd > heatingConditions =
            statePropagationFunction( orbital_element_conversions::convertKeplerianToCartesianElements(
                                          initialKeplerianState, planetaryGravitationalParameter ) ).second.second;

    // Separate time and heat rate and find maximum heat rate
    unsigned int i = 0;
    Eigen::VectorXd historyOfDynamicPressures;
    historyOfDynamicPressures.resize( heatingConditions.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = heatingConditions.begin( );
          mapIterator != heatingConditions.end( ); mapIterator++, i++ )
    {
        historyOfDynamicPressures[ i ] = mapIterator->second[ 0 ];
    }
    double dynamicPressure = historyOfDynamicPressures.maxCoeff( );

    // Return maximum value of dynamic pressure w.r.t. threshold value
    return dynamicPressure - minimumDynamicPressure;
}

//! Function to be used as input to the root-finder to determine the magnitude of the apoapsis maneuver.
double maneuverBisectionFunction( const double currentMagnitudeGuess, const Eigen::Vector6d& initialEstimatedCartesianState,
                                  const double targetPeriapsisRadius, const Eigen::Matrix3d& transformationFromLocalToInertialFrame,
                                  const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
                                  std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
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
    double predictedPeriapsisRadius = historyOfRadialDistances.minCoeff( );

    // Return predicted periapsis altitude w.r.t. expected value
    return predictedPeriapsisRadius - targetPeriapsisRadius;
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
                                         1.0 / currentParameterEstimate[ 1 ], { currentParameterEstimate[ 2 ],
                    currentParameterEstimate[ 3 ], currentParameterEstimate[ 4 ] } ) );

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

} // namespace system_models

} // namespace tudat
