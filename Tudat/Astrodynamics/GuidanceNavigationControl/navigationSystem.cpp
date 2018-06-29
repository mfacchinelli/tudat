#include <iostream>

#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to be used as input to the root-finder to determine the centroid of the acceleration curve.
double areaBisectionFunction( const double currentTrueAnomalyGuess, const std::vector< double >& estimatedTrueAnomaly,
                              const std::vector< double >& estimatedAerodynamicAccelerationMagnitude )
{
    // Find nearest lower index to true anomaly guess
    int nearestLowerIndex = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch( estimatedTrueAnomaly,
                                                                                            currentTrueAnomalyGuess );

    // Compute trapezoidal quadrature to integrate the area until and after the current guess
    double lowerSliceQuadratureResult = numerical_quadrature::performTrapezoidalQuadrature(
                utilities::sliceStlVector( estimatedTrueAnomaly, 0, nearestLowerIndex ),
                utilities::sliceStlVector( estimatedAerodynamicAccelerationMagnitude, 0, nearestLowerIndex ) );
    double upperSliceQuadratureResult = numerical_quadrature::performTrapezoidalQuadrature(
                utilities::sliceStlVector( estimatedTrueAnomaly, nearestLowerIndex + 1 ),
                utilities::sliceStlVector( estimatedAerodynamicAccelerationMagnitude, nearestLowerIndex + 1 ) );

    // Return difference in areas
    return upperSliceQuadratureResult - lowerSliceQuadratureResult;
}

//! Function to run the State Estimator (SE).
void NavigationSystem::runStateEstimator( const double previousTime, const Eigen::Vector4d& currentExternalMeasurementVector )
{
    // Save old true anomaly estimate
    double oldEstimatedTrueAnomaly = currentEstimatedKeplerianState_[ 5 ];

    // Update filter
    navigationFilter_->updateFilter( currentTime_, currentExternalMeasurementVector );

    // Extract estimated state and update navigation estimates
    Eigen::Vector16d updatedEstimatedState = navigationFilter->getCurrentStateEstimate( );
    setCurrentEstimatedCartesianState( updatedEstimatedState.segment( 0, 6 ) );
    currentEstimatedRotationalState_ = updatedEstimatedState.segment( 6, 4 );
    storeCurrentTimeAndStateEstimates( );

    // Check if new orbit and store new state estimate
    if ( ( oldEstimatedTrueAnomaly < 2.0 * mathematical_constants::PI ) &&
         ( currentEstimatedKeplerianState_[ 5 ] > 0 ) )
    {
        currentOrbitCounter_++;
    }


    // Store initial time and state and update body and acceleration maps
    storeCurrentTimeAndTranslationalStateEstimates( );
//    updateBodyAndAccelerationMaps( currentTime_ ); // done in onboard computer model
}

//! Function to run the Periapse Time Estimator (PTE).
void NavigationSystem::runPeriapseTimeEstimator( const std::map< double, Eigen::Vector3d >& mapOfEstimatedAerodynamicAcceleration )
{
    // Extract aerodynamic accelerations of when the spacecraft is below the atmospheric interface altitude
    double currentIterationTime;
    std::map< double, Eigen::Vector6d > mapOfEstimatedKeplerianStatesBelowAtmosphericInterface;
    std::vector< double > vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface;
    for ( std::map< double, std::pair< Eigen::VectorXd, Eigen::VectorXd > >::iterator stateIterator =
          currentOrbitHistoryOfEstimatedTranslationalStates_.begin( );
          stateIterator != currentOrbitHistoryOfEstimatedTranslationalStates_.end( ); stateIterator++ )
    {
        currentIterationTime = stateIterator->first;
        if ( stateIterator->second.first.segment( 0, 3 ).norm( ) <= atmosphericInterfaceRadius_ )
        {
            // Retireve time, state and acceleration of where the altitude is below the atmospheric interface
            mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ] =
                    stateIterator->second.second;
            vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.push_back(
                        mapOfEstimatedAerodynamicAcceleration[ currentIterationTime ].norm( ) );

            // Modify the true anomaly such that it is negative where it is above PI radians (before estimated periapsis)
            if ( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ][ 5 ] >= mathematical_constants::PI )
            {
                mapOfEstimatedKeplerianStatesBelowAtmosphericInterface[ currentIterationTime ][ 5 ] -= 2.0 * mathematical_constants::PI;
            }
        }
    }

    // Separate time and accelerations
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > pairOfEstimatedKeplerianStateBelowAtmosphericInterface =
            utilities::extractKeyAndValuesFromMap( mapOfEstimatedKeplerianStatesBelowAtmosphericInterface );
    Eigen::VectorXd timeBelowAtmosphericInterface = pairOfEstimatedKeplerianStateBelowAtmosphericInterface.first;
    Eigen::MatrixXd estimatedKeplerianStateBelowAtmosphericInterface = pairOfEstimatedKeplerianStateBelowAtmosphericInterface.second;
    Eigen::VectorXd estimatedTrueAnomalyBelowAtmosphericInterface = estimatedKeplerianStateBelowAtmosphericInterface.row( 5 ).transpose( );

    // Check that true anomaly and aerodynamic acceleration have the same length
    if ( estimatedTrueAnomalyBelowAtmosphericInterface.rows( ) !=
         vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface.size( ) )
    {
        throw std::runtime_error( "Error in periapse time estimator. The sizes of the true anomaly and aerodynamic accelerations "
                                  "vectors do not match." );
    }

    // Set root-finder boundaries as the first and last true anomaly elements
    areaBisectionRootFinder_->resetBoundaries(
                estimatedTrueAnomalyBelowAtmosphericInterface[ 0 ],
            estimatedTrueAnomalyBelowAtmosphericInterface[ estimatedTrueAnomalyBelowAtmosphericInterface.rows( ) - 1 ] );

    // Get intermediate variables
    Eigen::Vector6d initialEstimatedKeplerianState = estimatedKeplerianStateBelowAtmosphericInterface.col( 0 );

    // Set root-finder function as the area below the acceleration curve
    double estimatedErrorInTrueAnomaly = areaBisectionRootFinder_->execute(
                basic_mathematics::FunctionProxy(
                    boost::bind( &areaBisectionFunction, _1,
                                 utilities::convertEigenVectorToStlVector( estimatedTrueAnomalyBelowAtmosphericInterface ),
                                 vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface ) ) );
    std::cout << "Estimated Error in True Anomaly: " <<
                 unit_conversions::convertDegreesToRadians( estimatedErrorInTrueAnomaly ) << " deg" << std::endl;
    // note that this represents directly the error in estimated true anomaly, since the true anomaly of
    // periapsis is zero by definition

    // Find nearest lower index to error in true anomaly
    int estimatedPeriapsisIndex = basic_mathematics::computeNearestLeftNeighborUsingBinarySearch(
                estimatedErrorInTrueAnomaly, estimatedTrueAnomalyBelowAtmosphericInterface );

    // Compute estimated change in velocity (i.e., Delta V) due to aerodynamic acceleration
    double esimatedChangeInVelocity = - numerical_quadrature::performTrapezoidalQuadrature(
                utilities::convertEigenVectorToStlVector( timeBelowAtmosphericInterface ),
                vectorOfEstimatedAerodynamicAccelerationMagnitudeBelowAtmosphericInterface );
    std::cout << "Estimated Change in Velocity: " << esimatedChangeInVelocity << " m/s" << std::endl;

    // Compute estimated mean motion by using the semi-major axis at beginning of atmospheric phase
    double currentEstimatedMeanMotion = std::sqrt( planetaryGravitationalParameter_ /
                                                   std::pow( initialEstimatedKeplerianState[ 0 ], 3 ) );

    // Compute estimated change in semi-major axis due to change in velocity
    // This value is computed by assuming that the change in velocity occurs at pericenter and is given by
    // by an impulsive shot. Also here, the eccentricity at the beginning of the atmospheric phase is used
    double estimatedChangeInSemiMajorAxisDueToChangeInVelocity = 2.0 / currentEstimatedMeanMotion * std::sqrt(
                ( 1.0 + initialEstimatedKeplerianState[ 1 ] ) /
                ( 1.0 - initialEstimatedKeplerianState[ 1 ] ) ) * estimatedChangeInVelocity;
    std::cout << "Estimated Change in Semi-major Axis: " << estimatedChangeInSemiMajorAxis << " m" << std::endl;

    // Get estimated chagne in semi-major axis from estimated Keplerian state and find estimated error in semi-major axis
    double estimatedChangeInSemiMajorAxisFromKeplerianStateHistory = initialEstimatedKeplerianState[ 0 ] -
            estimatedKeplerianStateBelowAtmosphericInterface( 0, estimatedKeplerianStateBelowAtmosphericInterface.cols( ) );
    double estimatedErrorInSemiMajorAxis = estimatedChangeInSemiMajorAxisFromKeplerianStateHistory -
            estimatedChangeInSemiMajorAxisDueToChangeInVelocity;
    std::cout << "Estimated Error in Semi-major Axis: " << estimatedErrorInSemiMajorAxis << " m" << std::endl;

    // Combine errors to produce a vector of estimated error in Keplerian state
    Eigen::Vector6d estimatedErrorInKeplerianState = Eigen::Vector6d::Zero( );
    estimatedErrorInKeplerianState[ 0 ] = estimatedErrorInSemiMajorAxis;
    estimatedErrorInKeplerianState[ 5 ] = estimatedErrorInTrueAnomaly;

    // Correct latest estimated Keplerian state with new information from PTE
    double timeFromPeriapsisToCurrentPosition = timeBelowAtmosphericInterface[ timeBelowAtmosphericInterface.rows( ) ] -
            timeBelowAtmosphericInterface[ estimatedPeriapsisIndex ];
    Eigen::Vector6d updatedCurrentEstimatedKeplerianState = initialEstimatedKeplerianState;
    initialEstimatedKeplerianState[ 5 ] = 0.0; // set true anomaly to periapsis, to assume that no change in orbital elements
                                               // has occurred during atmospheric phase
    updatedCurrentEstimatedKeplerianState += stateTransitionMatrixFunction_( initialEstimatedKeplerianState ) *
            estimatedErrorInKeplerianState * timeFromPeriapsisToCurrentPosition;

    // Update navigation system state estimates
    setCurrentEstimatedKeplerianState( updatedCurrentEstimatedKeplerianState );
    historyOfEstimatedErrorsInKeplerianState_[ currentOrbitCounter_ ] = estimatedErrorInKeplerianState;

    // Update navigation filter state and covariance
    Eigen::Vector16d updatedCurrentEstimatedState = navigationFilter_->getCurrentStateEstimate( );
    updatedCurrentEstimatedState.segment( 0, 6 ) = currentEstimatedCartesianState_;
    navigationFilter_->modifyCurrentStateAndCovarianceEstimates( updatedCurrentEstimatedState,
                                                                 Eigen::Matrix16d::Identity( ) );
    // the covariance matrix is reset to the identity, since the new state is improved in accuracy
}

//! Function to run the Atmosphere Estimator (AE).
void NavigationSystem::runAtmosphereEstimator( const std::map< double, Eigen::Vector3d >& mapOfEstimatedAerodynamicAcceleration )
{

    // Update atmosphere settings of onboard body map
//    onboardBodyMap_.at( planetName_ )->setAtmosphereModel( );
}

//! Function to create navigation filter object for onboard state estimation.
void NavigationSystem::createNavigationFilter(
        const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardSystemModel,
        const boost::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) >& onboardMeasurementModel )
{
    // Create filter object
    navigationFilter_ = filters::createFilter( navigationFilterSettings_, &onboardSystemModel, &onboardMeasurementModel );

    // Set initial time
    currentTime_ = navigationFilter_->getInitialTime( );
    currentOrbitCounter_ = 0;

    // Retrieve navigation filter step size
    navigationRefreshStepSize_ = navigationFilter_->getIntegrationStepSize( );

    // Set initial state
    currentEstimatedCartesianState_ = navigationFilter_->getCurrentStateEstimate( );
    currentEstimatedKeplerianState_ =
            orbital_element_conversions::convertCartesianToKeplerianElements( currentEstimatedCartesianState_,
                                                                              planetaryGravitationalParameter_ );

    // Store initial time and state and update body and acceleration maps
    storeCurrentTimeAndTranslationalStateEstimates( );
    updateBodyAndAccelerationMaps( currentTime_ );
}

} // namespace navigation

} // namespace tudat
