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
double correctionFactorForCorridorBoundaries( const double altitudeGuess, const Eigen::Vector2d& linearLeastSquaresEstimate )
{
    return 1.0 / ( linearLeastSquaresEstimate[ 0 ] + linearLeastSquaresEstimate[ 1 ] * altitudeGuess );
}

//! Function to run the estimator for the specified corridor boundary.
double CorridorEstimator::estimateCorridorBoundary(
        const BoundaryType typeOfBoundaryToEstimate, const Eigen::Vector6d& initialEstimatedKeplerianState,
        const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
        std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
{
    // Set root-finder boundaries and compute altitude value
    double estimatedAltitudeBound;
    switch ( typeOfBoundaryToEstimate )
    {
    case lower_heating:
        // Set boundaries
        altitudeBisectionRootFinder_->resetBoundaries( boundariesForLowerAltitudeBasedOnHeating_.first,
                                                       boundariesForLowerAltitudeBasedOnHeating_.second );

        // Compute estimate via bisection root-finder
        estimatedAltitudeBound = altitudeBisectionRootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                        boost::bind( &CorridorEstimator::lowerAltitudeBisectionFunctionBasedOnHeatingConditions, this, _1,
                                     initialEstimatedKeplerianState, statePropagationFunction ) ) );
        break;
    case lower_lifetime:
        // Set boundaries
        altitudeBisectionRootFinder_->resetBoundaries( boundariesForLowerAltitudeBasedOnLifetime_.first,
                                                       boundariesForLowerAltitudeBasedOnLifetime_.second );

        // Compute estimate via bisection root-finder
        estimatedAltitudeBound = altitudeBisectionRootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                        boost::bind( &CorridorEstimator::lowerAltitudeBisectionFunctionBasedOnLifetimeCondition, this, _1,
                                     initialEstimatedKeplerianState, statePropagationFunction ) ) );
        break;
    case upper:
        // Set boundaries
        altitudeBisectionRootFinder_->resetBoundaries( boundariesForUpperAltitude_.first, boundariesForUpperAltitude_.second );

        // Compute estimate via bisection root-finder
        estimatedAltitudeBound = altitudeBisectionRootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                        boost::bind( &CorridorEstimator::upperAltitudeBisectionFunction, this, _1,
                                     initialEstimatedKeplerianState, statePropagationFunction ) ) );
        break;
    }

    // Give output
    return estimatedAltitudeBound;
}

//! Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
double CorridorEstimator::lowerAltitudeBisectionFunctionBasedOnHeatingConditions(
        const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
        const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
        std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
{
    // Propagate orbit to new condition and retrieve heating conditions
    std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationResult =
            propagateStateWithAltitudeGuess( currentAltitudeGuess, initialEstimatedKeplerianState, statePropagationFunction );
    std::map< double, Eigen::VectorXd > trajectory = propagationResult.first;
    std::map< double, Eigen::VectorXd > heatingConditions = propagationResult.second;

    // Separate time and heat rate and find maximum heat rate
    unsigned int i = 0;
    std::vector< double > historyOfTimes;
    Eigen::VectorXd historyOfHeatRates;
    historyOfHeatRates.resize( heatingConditions.size( ) );
    Eigen::VectorXd historyOfAltitudes;
    historyOfAltitudes.resize( trajectory.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = heatingConditions.begin( );
          mapIterator != heatingConditions.end( ); mapIterator++, i++ )
    {
        historyOfTimes.push_back( mapIterator->first );
        historyOfHeatRates[ i ] = mapIterator->second[ 1 ];
        historyOfAltitudes[ i ] = trajectory[ mapIterator->first ].segment( 0, 3 ).norm( ) - planetaryRadius_;
    }
    double heatRate = historyOfHeatRates.maxCoeff( );

    // Compute heat load by integrating heat flux
    double heatLoad = numerical_quadrature::performTrapezoidalQuadrature(
                historyOfTimes, utilities::convertEigenVectorToStlVector( historyOfHeatRates ) );

    // Compute offsets w.r.t. maximum allowed heat rate and heat load
    double offsetInHeatRate = heatRate - maximumHeatRate_;
    double offsetInHeatLoad = heatLoad - maximumHeatLoad_;

    // Extract and store actual periapsis
    historyOfPeriapsisInformation_.push_back( std::make_pair( currentAltitudeGuess, historyOfAltitudes.minCoeff( ) ) );

    // Return value to to indicate closeness to limiting value
    return ( offsetInHeatRate > offsetInHeatLoad ) ? offsetInHeatRate : offsetInHeatLoad;
}

//! Function to be used as input to the root-finder to determine the lower altitude bound for the periapsis corridor.
double CorridorEstimator::lowerAltitudeBisectionFunctionBasedOnLifetimeCondition(
        const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
        const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
        std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
{
    // Propagate orbit to new condition and retrieve heating conditions
    std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationResult =
            propagateStateWithAltitudeGuess( currentAltitudeGuess, initialEstimatedKeplerianState, statePropagationFunction );
    std::map< double, Eigen::VectorXd > trajectory = propagationResult.first;

    // Extract history of altitudes
    unsigned int i = 0;
    Eigen::VectorXd historyOfAltitudes;
    historyOfAltitudes.resize( trajectory.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = trajectory.begin( );
          mapIterator != trajectory.end( ); mapIterator++, i++ )
    {
        historyOfAltitudes[ i ] = mapIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius_;
    }

    // Extract and store actual periapsis
    historyOfPeriapsisInformation_.push_back( std::make_pair( currentAltitudeGuess, historyOfAltitudes.minCoeff( ) ) );

    // Give propagation a score based on lifetime
    double actualLifetime = ( trajectory.rbegin( )->first - trajectory.begin( )->first ) / physical_constants::JULIAN_DAY;

    // Return value to to indicate closeness to limiting value
    return actualLifetime - minimumLifetime_;
}

//! Function to be used as input to the root-finder to determine the upper altitude bound for the periapsis corridor.
double CorridorEstimator::upperAltitudeBisectionFunction(
        const double currentAltitudeGuess, const Eigen::Vector6d& initialEstimatedKeplerianState,
        const boost::function< std::pair< bool, std::pair< std::map< double, Eigen::VectorXd >,
        std::map< double, Eigen::VectorXd > > >( const Eigen::Vector6d& ) >& statePropagationFunction )
{
    // Propagate orbit to new condition and retrieve heating conditions
    std::pair< std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd > > propagationResult =
            propagateStateWithAltitudeGuess( currentAltitudeGuess, initialEstimatedKeplerianState, statePropagationFunction );
    std::map< double, Eigen::VectorXd > trajectory = propagationResult.first;
    std::map< double, Eigen::VectorXd > heatingConditions = propagationResult.second;

    // Separate time and heat rate and find maximum heat rate
    unsigned int i = 0;
    Eigen::VectorXd historyOfDynamicPressures;
    historyOfDynamicPressures.resize( heatingConditions.size( ) );
    Eigen::VectorXd historyOfAltitudes;
    historyOfAltitudes.resize( trajectory.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = heatingConditions.begin( );
          mapIterator != heatingConditions.end( ); mapIterator++, i++ )
    {
        historyOfDynamicPressures[ i ] = mapIterator->second[ 0 ];
        historyOfAltitudes[ i ] = trajectory[ mapIterator->first ].segment( 0, 3 ).norm( ) - planetaryRadius_;
    }
    double dynamicPressure = historyOfDynamicPressures.maxCoeff( );

    // Extract and store actual periapsis
    historyOfPeriapsisInformation_.push_back( std::make_pair( currentAltitudeGuess, historyOfAltitudes.minCoeff( ) ) );

    // Return maximum value of dynamic pressure w.r.t. threshold value
    return dynamicPressure - minimumDynamicPressure_;
}

} // namespace system_models

} // namespace tudat
