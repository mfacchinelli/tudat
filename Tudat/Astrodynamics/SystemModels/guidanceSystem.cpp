#include "Tudat/Astrodynamics/SystemModels/guidanceSystem.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/SystemModels/extraFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"

namespace tudat
{

namespace system_models
{

//! Function to run corridor estimator (CE).
void GuidanceSystem::runCorridorEstimator( const double currentTime,
                                           const Eigen::Vector6d& currentEstimatedCartesianState,
                                           const Eigen::Vector6d& currentEstimatedKeplerianState,
                                           const double planetaryRadius,
                                           const double planetaryGravitationalParameter )
{
    // Inform user
    std::cout << std::endl << "Estimating Periapsis Corridor." << std::endl;

    // Loop until convergence is reached
    bool corridorEstimationComplete = true;
    do
    {
        // Run corridor estimator with heating conditions as lower limit
        estimateCorridorBoundaries( currentTime, currentEstimatedCartesianState, currentEstimatedKeplerianState,
                                    planetaryRadius, planetaryGravitationalParameter );

        // If in walk-out phase, check that lifetime is above minimum value
        if ( currentOrbitAerobrakingPhase_ == walk_out_phase )
        {
            // Inform user
            std::cout << "Walk-out phase. Checking lifetime constraints." << std::endl;

            // Check that lifetime threshold is met
            corridorEstimationComplete = lifetimeReducedStatePropagationFunction_( currentEstimatedCartesianState ).first;
            if ( !corridorEstimationComplete )
            {
                // Inform user
                std::cout << "Lifetime requirement not met. Re-estimating corridor boundaries." << std::endl;

                // Run corridor estimator with lifetime as lower limit
                estimateCorridorBoundaries( currentTime, currentEstimatedCartesianState, currentEstimatedKeplerianState,
                                            planetaryRadius, planetaryGravitationalParameter, false );
                corridorEstimationComplete = true;
            }
        }
    }
    while ( !corridorEstimationComplete );
}

//! Function to run maneuver estimator (ME).
void GuidanceSystem::runManeuverEstimator( const Eigen::Vector6d& currentEstimatedCartesianState,
                                           const Eigen::Vector6d& currentEstimatedKeplerianState,
                                           const double currentEstimatedMeanMotion,
                                           const double planetaryRadius,
                                           const bool improveEstimateWithBisection )
{
    // Inform user
    std::cout << std::endl << "Estimating Apoapsis Maneuver." << std::endl;

    // Set apoapsis maneuver vector to zero
    scheduledApsoapsisManeuver_.setZero( );

    // Compute predicted periapsis radius
    double predictedPeriapsisAltitude = computeCurrentFirstOrderEstimatedPeriapsisRadius( currentEstimatedKeplerianState ) - planetaryRadius;
    double differenceInPeriapsisAltitude = std::get< 2 >( periapsisTargetingInformation_ ) - predictedPeriapsisAltitude;

    // Compute estimated maneuver in y-direction of local orbit frame
    double preliminaryApoapsisManeuverMagnitude = 0.25 * currentEstimatedMeanMotion * differenceInPeriapsisAltitude * std::sqrt(
                ( 1.0 + currentEstimatedKeplerianState[ 1 ] ) / ( 1.0 - currentEstimatedKeplerianState[ 1 ] ) );
    std::cout << "Preliminary magnitude: " << preliminaryApoapsisManeuverMagnitude << " m/s" << std::endl;

    // Compute transformation from local to inertial frame
    Eigen::Matrix3d transformationFromLocalToInertialFrame =
            computeCurrentRotationFromLocalToInertialFrame( currentEstimatedCartesianState );

    // Set root-finder boundaries for the maneuver estimation
    maneuverBisectionRootFinder_->resetBoundaries( 2.0 * preliminaryApoapsisManeuverMagnitude, 0.5 * preliminaryApoapsisManeuverMagnitude );

    // Improve the estimate if magnitude is large enough
    double estimatedApoapsisManeuverMagnitude;
    if ( improveEstimateWithBisection && ( std::fabs( preliminaryApoapsisManeuverMagnitude ) > 0.15 ) ) // 0.15 N is an empirical value
    {
        // Try using root-finder to improve estimate
        try
        {
            // Set root-finder function as the heat rate and heat load calculator
            estimatedApoapsisManeuverMagnitude = maneuverBisectionRootFinder_->execute(
                        boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                            boost::bind( &maneuverBisectionFunction, _1, currentEstimatedCartesianState,
                                         std::get< 2 >( periapsisTargetingInformation_ ) + planetaryRadius,
                                         transformationFromLocalToInertialFrame, periodReducedStatePropagationFunction_ ) ) );
            std::cout << "Improved estimate: " << estimatedApoapsisManeuverMagnitude << " m/s" << std::endl
                      << "Ratio: " << estimatedApoapsisManeuverMagnitude / preliminaryApoapsisManeuverMagnitude << std::endl;
        }
        catch ( std::runtime_error& caughtException )
        {
            // Inform user on error
            std::cerr << "Error while computing improved estimate for apoapsis maneuver. Caught this exception during root-finder "
                         "operation: " << caughtException.what( ) << std::endl
                      << "The preliminary magnitude will be used to carry out the maneuver." << std::endl;

            // Take preliminary magnitude as maneuver magnitude
            estimatedApoapsisManeuverMagnitude = preliminaryApoapsisManeuverMagnitude;
        }
    }
    else
    {
        estimatedApoapsisManeuverMagnitude = preliminaryApoapsisManeuverMagnitude;
    }

    // Add magnitude to maneuver vector and save to hisotry
    scheduledApsoapsisManeuver_[ 1 ] = estimatedApoapsisManeuverMagnitude;
    historyOfApoapsisManeuverMagnitudes_[ currentOrbitCounter_ ] = estimatedApoapsisManeuverMagnitude;

    // Convert manveuver (i.e., Delta V vector) to inertial frame
    scheduledApsoapsisManeuver_ = transformationFromLocalToInertialFrame * scheduledApsoapsisManeuver_;
    std::cout << "Scheduled maneuver: " << scheduledApsoapsisManeuver_.transpose( ) << std::endl;
}

//! Function to run corridor estimator with nominal conditions.
void GuidanceSystem::estimateCorridorBoundaries( const double currentTime,
                                                 const Eigen::Vector6d& currentEstimatedCartesianState,
                                                 const Eigen::Vector6d& currentEstimatedKeplerianState,
                                                 const double planetaryRadius,
                                                 const double planetaryGravitationalParameter,
                                                 const bool useHeatAsLowerBoundaryThreshold )
{
    // Create propagation termination settings based on period and lifetime to be used throughout the function
    double periodTerminationTime = currentTime + 2.0 / 3.0 *
            basic_astrodynamics::computeKeplerOrbitalPeriod( currentEstimatedKeplerianState[ 0 ], planetaryGravitationalParameter );
    boost::shared_ptr< propagators::PropagationTerminationSettings > periodTerminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( periodTerminationTime );
    double lifetimeTerminationTime = currentTime + 1.5 * minimumAllowedLifetime_ * physical_constants::JULIAN_DAY;
    boost::shared_ptr< propagators::PropagationTerminationSettings > lifetimeTerminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( lifetimeTerminationTime );

    // Create reduced state propagation function where termination settings are already set
    periodReducedStatePropagationFunction_ = boost::bind( statePropagationFunction_, periodTerminationSettings, _1 );
    lifetimeReducedStatePropagationFunction_ = boost::bind( statePropagationFunction_, lifetimeTerminationSettings, _1 );

    // Propagate state for two thirds of the orbit
    std::map< double, Eigen::VectorXd > unaffectedPropagatedState =
            periodReducedStatePropagationFunction_( currentEstimatedCartesianState ).second.first;

    // Retrieve periapsis altitude
    unsigned int i = 0;
    Eigen::VectorXd historyOfAltitudes;
    historyOfAltitudes.resize( unaffectedPropagatedState.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = unaffectedPropagatedState.begin( );
          mapIterator != unaffectedPropagatedState.end( ); mapIterator++, i++ )
    {
        historyOfAltitudes[ i ] = mapIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
    }
    double predictedPeriapsisAltitude = historyOfAltitudes.minCoeff( );
    std::cout << "Predicted periapsis altitude: " << predictedPeriapsisAltitude / 1.0e3 << " km" << std::endl;

    // Compute lower bound of corridor based on heat or lifetime
    double estimatedLowerAltitudeBound;
    if ( useHeatAsLowerBoundaryThreshold )
    {
        // Set root-finder boundaries for the lower altitude limits
        altitudeBisectionRootFinder_->resetBoundaries( 90.0e3, 120.0e3 );

        // Set root-finder function as the heat rate and heat load calculator
        // Add scaling if spacecraft is in walk-in or walk-out phases
        estimatedLowerAltitudeBound = periapsisAltitudeScaling_ * altitudeBisectionRootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                        boost::bind( &lowerAltitudeBisectionFunctionBasedOnHeatingConditions, _1, currentEstimatedKeplerianState,
                                     planetaryRadius, planetaryGravitationalParameter, maximumAllowedHeatRate_, maximumAllowedHeatLoad_,
                                     periodReducedStatePropagationFunction_ ) ) );
    }
    else
    {
        // Set root-finder boundaries for the lower altitude limits
        altitudeBisectionRootFinder_->resetBoundaries( 90.0e3, 150.0e3 );

        // Set root-finder function as the lifetime calculator
        estimatedLowerAltitudeBound = altitudeBisectionRootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                        boost::bind( &lowerAltitudeBisectionFunctionBasedOnLifetimeCondition, _1, currentEstimatedKeplerianState,
                                     planetaryRadius, planetaryGravitationalParameter, minimumAllowedLifetime_,
                                     lifetimeReducedStatePropagationFunction_ ) ) );
    }
    std::cout << "Lower boundary: " << estimatedLowerAltitudeBound / 1.0e3 << " km" << std::endl;

    // Set root-finder boundaries for the upper altitude limits
    altitudeBisectionRootFinder_->resetBoundaries( 100.0e3, 150.0e3 );

    // Set root-finder function as the dynamic pressure calculator
    // Add scaling if spacecraft is in walk-in or walk-out phases
    double estimatedUpperAltitudeBound = periapsisAltitudeScaling_ * altitudeBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                    boost::bind( &upperAltitudeBisectionFunction, _1, currentEstimatedKeplerianState, planetaryRadius,
                                 planetaryGravitationalParameter, minimumAllowedDynamicPressure_, periodReducedStatePropagationFunction_ ) ) );
    std::cout << "Upper boundary: " << estimatedUpperAltitudeBound / 1.0e3 << " km" << std::endl;

    // Check that lower altitude bound is indeed lower
    if ( estimatedLowerAltitudeBound > estimatedUpperAltitudeBound )
    {
        // Inform user of altitude conflict
        std::cerr << "Warning in periapsis corridor estimator. The lower altitude bound is larger than the higher altitude bound. "
                     "The upper bound will be defined as the lower bound plus 15 km." << std::endl;

        // Replace upper altitude with periapsis estimate
        estimatedUpperAltitudeBound = estimatedLowerAltitudeBound + 15.0e3;
    }

    // Check whether predicted altitude is within bounds
    bool isAltitudeWithinBounds = ( ( predictedPeriapsisAltitude > estimatedLowerAltitudeBound ) &&
                                    ( predictedPeriapsisAltitude < estimatedUpperAltitudeBound ) );

    // Compute target altitude
    // If the predicted periapsis altitude falls within the bounds, then the target altitude is simply the
    // predicted one, otherwise it is the mid-point of the corridor
    double estimatedTargetPeriapsisAltitude = isAltitudeWithinBounds ? predictedPeriapsisAltitude :
                                                                       0.5 * ( estimatedLowerAltitudeBound +
                                                                               estimatedUpperAltitudeBound );
    std::cout << "Target periapsis altitude: " << estimatedTargetPeriapsisAltitude / 1.0e3 << " km" << std::endl;

    // Save periapsis corridor altitudes to history
    historyOfEstimatedPeriapsisCorridorBoundaries_[ currentOrbitCounter_ ] =
            std::make_pair( estimatedLowerAltitudeBound, estimatedUpperAltitudeBound );

    // Store corridor information to pair
    periapsisTargetingInformation_ = std::make_tuple( isAltitudeWithinBounds, predictedPeriapsisAltitude,
                                                      estimatedTargetPeriapsisAltitude );
}

} // namespace system_models

} // namespace tudat
