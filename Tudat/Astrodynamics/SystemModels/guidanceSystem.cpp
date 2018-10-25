#include "Tudat/Astrodynamics/SystemModels/guidanceSystem.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
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
    // Run corridor estimator if aerobraking is not complete
    switch ( currentOrbitAerobrakingPhase_ )
    {
    case aerobraking_complete:
    {
        // Inform user
        std::cout << std::endl << "Aerobraking Complete." << std::endl;

        // Save periapsis corridor altitudes to history
        historyOfEstimatedPeriapsisCorridorBoundaries_[ currentOrbitCounter_ ] = std::make_pair( TUDAT_NAN, TUDAT_NAN );

        // Store corridor information to pair
        periapsisTargetingInformation_ = std::make_tuple( true, targetPeriapsisAltitude_, targetPeriapsisAltitude_ );
        break;
    }
    case termination_phase:
    {
        // Inform user
        std::cout << std::endl << "Raising Periapsis To Target Value." << std::endl;

        // Create propagation termination settings based on period and lifetime to be used throughout the function
        double periodTerminationTime = currentTime + 2.0 / 3.0 *
                basic_astrodynamics::computeKeplerOrbitalPeriod( currentEstimatedKeplerianState[ 0 ], planetaryGravitationalParameter );
        boost::shared_ptr< propagators::PropagationTerminationSettings > periodTerminationSettings =
                boost::make_shared< propagators::PropagationTimeTerminationSettings >( periodTerminationTime );

        // Propagate state for two thirds of the orbit
        std::map< double, Eigen::VectorXd > nominalPropagatedState =
                statePropagationFunction_( periodTerminationSettings, currentEstimatedCartesianState ).second.first;

        // Retrieve periapsis altitude
        unsigned int i = 0;
        Eigen::VectorXd historyOfAltitudes;
        historyOfAltitudes.resize( nominalPropagatedState.size( ) );
        for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = nominalPropagatedState.begin( );
              mapIterator != nominalPropagatedState.end( ); mapIterator++, i++ )
        {
            historyOfAltitudes[ i ] = mapIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
        }
        double predictedPeriapsisAltitude = historyOfAltitudes.minCoeff( );
        std::cout << "Predicted periapsis altitude: " << predictedPeriapsisAltitude / 1.0e3 << " km" << std::endl
                  << "Target periapsis altitude: " << targetPeriapsisAltitude_ / 1.0e3 << " km" << std::endl;

        // Save periapsis corridor altitudes to history
        historyOfEstimatedPeriapsisCorridorBoundaries_[ currentOrbitCounter_ ] = std::make_pair( TUDAT_NAN, TUDAT_NAN );

        // Store corridor information to pair
        periapsisTargetingInformation_ = std::make_tuple( false, predictedPeriapsisAltitude, targetPeriapsisAltitude_ );
        break;
    }
    default:
    {
        // Inform user
        std::cout << std::endl << "Estimating Periapsis Corridor." << std::endl;

        // Run corridor estimator with heating conditions as lower limit
        estimateCorridorBoundaries( currentTime, currentEstimatedCartesianState, currentEstimatedKeplerianState,
                                    planetaryRadius, planetaryGravitationalParameter );

        // If in walk-out phase, check that lifetime is above minimum value
        if ( currentOrbitAerobrakingPhase_ == walk_out_phase )
        {
            // Inform user
            std::cout << "Walk-out phase. Checking lifetime constraints." << std::endl;

            // Check that lifetime threshold is met
            if ( !lifetimeReducedStatePropagationFunction_( currentEstimatedCartesianState ).first )
            {
                // Inform user
                std::cout << "Lifetime requirement not met. Re-estimating corridor boundaries." << std::endl;

                // Run corridor estimator with lifetime as lower limit
                estimateCorridorBoundaries( currentTime, currentEstimatedCartesianState, currentEstimatedKeplerianState,
                                            planetaryRadius, planetaryGravitationalParameter, false );
            }
            else
            {
                // Inform user
                std::cout << "Lifetime requirement met." << std::endl;
            }
        }
        break;
    }
    }
}

//! Function to run apoapsis maneuver estimator (ME).
void GuidanceSystem::runApoapsisManeuverEstimator( const Eigen::Vector6d& currentEstimatedCartesianState,
                                                   const Eigen::Vector6d& currentEstimatedKeplerianState,
                                                   const double currentEstimatedMeanMotion,
                                                   const double planetaryRadius,
                                                   const bool improveEstimateWithBisection )
{
    // Inform user
    std::cout << std::endl << "Estimating Apoapsis Maneuver." << std::endl;

    // Set apoapsis maneuver vector to zero
    scheduledApsisManeuver_.setZero( );

    // Compute difference in periapsis radius
    double differenceInPeriapsisAltitude = std::get< 2 >( periapsisTargetingInformation_ ) - std::get< 1 >( periapsisTargetingInformation_ );

    // Compute estimated maneuver in y-direction of local orbit frame
    double preliminaryApoapsisManeuverMagnitude = 0.25 * currentEstimatedMeanMotion * differenceInPeriapsisAltitude * std::sqrt(
                ( 1.0 + currentEstimatedKeplerianState[ 1 ] ) / ( 1.0 - currentEstimatedKeplerianState[ 1 ] ) );
    std::cout << "Preliminary magnitude: " << preliminaryApoapsisManeuverMagnitude << " m/s" << std::endl;

    // Compute transformation from local to inertial frame
    Eigen::Matrix3d transformationFromLocalToInertialFrame =
            computeCurrentRotationFromLocalToInertialFrame( currentEstimatedCartesianState );

    // Improve the estimate if magnitude is large enough
    double estimatedApoapsisManeuverMagnitude;
    if ( improveEstimateWithBisection &&
         ( std::fabs( preliminaryApoapsisManeuverMagnitude ) > 0.15 ) && // 0.15 N is an empirical value
         ( currentOrbitAerobrakingPhase_ != termination_phase ) )
    {
        // Try using root-finder to improve estimate
        try
        {
            // Set root-finder boundaries for the maneuver estimation
            maneuverBisectionRootFinder_->resetBoundaries( 2.0 * preliminaryApoapsisManeuverMagnitude,
                                                           0.5 * preliminaryApoapsisManeuverMagnitude );

            // Set root-finder function as the heat rate and heat load calculator
            estimatedApoapsisManeuverMagnitude = maneuverBisectionRootFinder_->execute(
                        boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                            boost::bind( &maneuverBisectionFunction, _1, currentEstimatedCartesianState,
                                         std::get< 2 >( periapsisTargetingInformation_ ) + planetaryRadius,
                                         transformationFromLocalToInertialFrame, periodReducedStatePropagationFunction_, true ) ) );
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
    scheduledApsisManeuver_[ 1 ] = estimatedApoapsisManeuverMagnitude;
    historyOfApsisManeuverMagnitudes_[ currentOrbitCounter_ ] = estimatedApoapsisManeuverMagnitude;

    // Convert manveuver (i.e., Delta V vector) to inertial frame
    scheduledApsisManeuver_ = transformationFromLocalToInertialFrame * scheduledApsisManeuver_;
    std::cout << "Scheduled maneuver: " << scheduledApsisManeuver_.transpose( ) << std::endl;
}

//! Function to run periapsis maneuver estimator (ME).
void GuidanceSystem::runPeriapsisManeuverEstimator( const double currentTime,
                                                    const Eigen::Vector6d& currentEstimatedCartesianState,
                                                    const Eigen::Vector6d& currentEstimatedKeplerianState,
                                                    const double currentEstimatedMeanMotion,
                                                    const double planetaryRadius,
                                                    const double planetaryGravitationalParameter )
{
    // Inform user
    std::cout << std::endl << "Estimating Periapsis Maneuver." << std::endl;

    // Create propagation termination settings based on period and lifetime to be used throughout the function
    double periodTerminationTime = currentTime + 2.0 / 3.0 *
            basic_astrodynamics::computeKeplerOrbitalPeriod( currentEstimatedKeplerianState[ 0 ], planetaryGravitationalParameter );
    boost::shared_ptr< propagators::PropagationTerminationSettings > periodTerminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( periodTerminationTime );

    // Create reduced state propagation functions where termination settings are already set
    periodReducedStatePropagationFunction_ = boost::bind( statePropagationFunction_, periodTerminationSettings, _1 );

    // Propagate state for two thirds of the orbit
    std::map< double, Eigen::VectorXd > nominalPropagatedState =
            periodReducedStatePropagationFunction_( currentEstimatedCartesianState ).second.first;

    // Retrieve periapsis altitude
    unsigned int i = 0;
    Eigen::VectorXd historyOfAltitudes;
    historyOfAltitudes.resize( nominalPropagatedState.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = nominalPropagatedState.begin( );
          mapIterator != nominalPropagatedState.end( ); mapIterator++, i++ )
    {
        historyOfAltitudes[ i ] = mapIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
    }
    double predictedApoapsisAltitude = historyOfAltitudes.maxCoeff( );
    std::cout << "Predicted apoapsis altitude: " << predictedApoapsisAltitude / 1.0e3 << " km" << std::endl
              << "Target apoapsis altitude: " << targetApoapsisAltitude_ / 1.0e3 << " km" << std::endl;

    // Set apoapsis maneuver vector to zero
    scheduledApsisManeuver_.setZero( );

    // Compute difference in periapsis radius
    double differenceInApoapsisAltitude = targetApoapsisAltitude_ - predictedApoapsisAltitude;

    // Compute estimated maneuver in y-direction of local orbit frame
    double preliminaryPeriapsisManeuverMagnitude = 0.25 * currentEstimatedMeanMotion * differenceInApoapsisAltitude * std::sqrt(
                ( 1.0 - currentEstimatedKeplerianState[ 1 ] ) / ( 1.0 + currentEstimatedKeplerianState[ 1 ] ) );
    std::cout << "Preliminary magnitude: " << preliminaryPeriapsisManeuverMagnitude << " m/s" << std::endl;

    // Compute transformation from local to inertial frame
    Eigen::Matrix3d transformationFromLocalToInertialFrame =
            computeCurrentRotationFromLocalToInertialFrame( currentEstimatedCartesianState );

    // Try using root-finder to improve estimate
    double estimatedPeriapsisManeuverMagnitude;
    try
    {
        // Set root-finder boundaries for the maneuver estimation
        maneuverBisectionRootFinder_->resetBoundaries( 2.0 * preliminaryPeriapsisManeuverMagnitude,
                                                       0.5 * preliminaryPeriapsisManeuverMagnitude );

        // Set root-finder function as the heat rate and heat load calculator
        estimatedPeriapsisManeuverMagnitude = maneuverBisectionRootFinder_->execute(
                    boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                        boost::bind( &maneuverBisectionFunction, _1, currentEstimatedCartesianState, targetApoapsisAltitude_ + planetaryRadius,
                                     transformationFromLocalToInertialFrame, periodReducedStatePropagationFunction_, false ) ) );
        std::cout << "Improved estimate: " << estimatedPeriapsisManeuverMagnitude << " m/s" << std::endl
                  << "Ratio: " << estimatedPeriapsisManeuverMagnitude / preliminaryPeriapsisManeuverMagnitude << std::endl;
    }
    catch ( std::runtime_error& caughtException )
    {
        // Inform user on error
        std::cerr << "Error while computing improved estimate for periapsis maneuver. Caught this exception during root-finder "
                     "operation: " << caughtException.what( ) << std::endl
                  << "The preliminary magnitude will be used to carry out the maneuver." << std::endl;

        // Take preliminary magnitude as maneuver magnitude
        estimatedPeriapsisManeuverMagnitude = preliminaryPeriapsisManeuverMagnitude;
    }

    // Add magnitude to maneuver vector and save to hisotry
    scheduledApsisManeuver_[ 1 ] = estimatedPeriapsisManeuverMagnitude;
    historyOfApsisManeuverMagnitudes_[ currentOrbitCounter_ ] = estimatedPeriapsisManeuverMagnitude;

    // Convert manveuver (i.e., Delta V vector) to inertial frame
    scheduledApsisManeuver_ = transformationFromLocalToInertialFrame * scheduledApsisManeuver_;
    std::cout << "Scheduled maneuver: " << scheduledApsisManeuver_.transpose( ) << std::endl;
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

    // Create reduced state propagation functions where termination settings are already set
    periodReducedStatePropagationFunction_ = boost::bind( statePropagationFunction_, periodTerminationSettings, _1 );
    lifetimeReducedStatePropagationFunction_ = boost::bind( statePropagationFunction_, lifetimeTerminationSettings, _1 );

    // Propagate state for two thirds of the orbit
    std::map< double, Eigen::VectorXd > nominalPropagatedState =
            periodReducedStatePropagationFunction_( currentEstimatedCartesianState ).second.first;

    // Retrieve periapsis altitude
    unsigned int i = 0;
    Eigen::VectorXd historyOfAltitudes;
    historyOfAltitudes.resize( nominalPropagatedState.size( ) );
    for ( std::map< double, Eigen::VectorXd >::const_iterator mapIterator = nominalPropagatedState.begin( );
          mapIterator != nominalPropagatedState.end( ); mapIterator++, i++ )
    {
        historyOfAltitudes[ i ] = mapIterator->second.segment( 0, 3 ).norm( ) - planetaryRadius;
    }
    double predictedPeriapsisAltitude = historyOfAltitudes.minCoeff( );
    std::cout << "Predicted periapsis altitude: " << predictedPeriapsisAltitude / 1.0e3 << " km" << std::endl;

    // Compute lower bound of corridor based on heat or lifetime
    double estimatedLowerAltitudeBound;
    if ( useHeatAsLowerBoundaryThreshold )
    {
        // Compute lower altitude bound based on heat rate and heat load
        // Add scaling if spacecraft is in walk-in phase
        estimatedLowerAltitudeBound = periapsisAltitudeScaling_ * corridorEstimator_->estimateCorridorBoundary(
                    CorridorEstimator::lower_heating, currentEstimatedKeplerianState, periodReducedStatePropagationFunction_ );
    }
    else
    {
        // Compute lower altitude bound based on lifetime
        // Add scaling if spacecraft is in walk-in phase
        estimatedLowerAltitudeBound = periapsisAltitudeScaling_ * corridorEstimator_->estimateCorridorBoundary(
                    CorridorEstimator::lower_lifetime, currentEstimatedKeplerianState, lifetimeReducedStatePropagationFunction_ );
    }

    // Compute upper bound of corridor based on dynamic pressure, or load from cache
    double estimatedUpperAltitudeBound;
    if ( useHeatAsLowerBoundaryThreshold )
    {
        // Compute upper altitude bound based on dynamic pressure
        // Add scaling if spacecraft is in walk-in phase
        estimatedUpperAltitudeBound = periapsisAltitudeScaling_ * corridorEstimator_->estimateCorridorBoundary(
                    CorridorEstimator::upper, currentEstimatedKeplerianState, periodReducedStatePropagationFunction_ );
    }
    else
    {
        // Use value from previous computation (the corridor estimator with heating conditions is always run first)
        estimatedUpperAltitudeBound = historyOfEstimatedPeriapsisCorridorBoundaries_[ currentOrbitCounter_ ].second;
    }

    // Correct result for difference between Keplerian assumption and reality
    altitudeCorrectionFunction_ = corridorEstimator_->getPeriapsisAltitudeCorrectionFunction( );
    estimatedLowerAltitudeBound *= altitudeCorrectionFunction_( estimatedLowerAltitudeBound );
    std::cout << "Lower alt corr: " << altitudeCorrectionFunction_( estimatedLowerAltitudeBound ) << std::endl;
    if ( useHeatAsLowerBoundaryThreshold ) // do not correct twice (if loaded from cache, correction has already been applied)
    {
        estimatedUpperAltitudeBound *= altitudeCorrectionFunction_( estimatedUpperAltitudeBound );
        std::cout << "Upper alt corr: " << altitudeCorrectionFunction_( estimatedUpperAltitudeBound ) << std::endl;
    }
    std::cout << "Lower boundary: " << estimatedLowerAltitudeBound / 1.0e3 << " km" << std::endl
              << "Upper boundary: " << estimatedUpperAltitudeBound / 1.0e3 << " km" << std::endl;

    // Check that lower altitude bound is indeed lower
    if ( estimatedLowerAltitudeBound > estimatedUpperAltitudeBound )
    {
        // Inform user of altitude conflict
        std::cerr << "Warning in periapsis corridor estimator. The lower altitude bound is larger than the higher altitude bound. "
                     "The upper bound will be defined as the lower bound plus 2.5 km." << std::endl;

        // Replace upper altitude with periapsis estimate
        estimatedUpperAltitudeBound = estimatedLowerAltitudeBound + 2.5e3;
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
