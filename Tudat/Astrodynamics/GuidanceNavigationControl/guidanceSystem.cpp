#include "Tudat/Astrodynamics/GuidanceNavigationControl/guidanceSystem.h"

#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/extraFunctions.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to run corridor estimator (CE).
void GuidanceSystem::runCorridorEstimator( const Eigen::Vector6d& currentEstimatedKeplerianState,
                                           const double planetaryRadius,
                                           const double planetaryGravitationalParameter,
                                           const double atmosphericInterfaceRadius,
                                           const boost::function< double( double ) >& densityFunction )
{
    // Inform user
    std::cout << "Estimating Periapsis Corridor." << std::endl;

    // Compute estimated apoapsis radius
    double estimatedApoapsisRadius = computeCurrentEstimatedApoapsisRadius( currentEstimatedKeplerianState );

    // Compute predicted periapsis altitude
    double predictedPeriapsisAltitude = computeCurrentEstimatedPeriapsisRadius( currentEstimatedKeplerianState ) - planetaryRadius;
    std::cout << "Predicted periapsis altitude: " << predictedPeriapsisAltitude / 1.0e3 << " km" << std::endl;

    // Set root-finder boundaries for the lower altitude limits
    altitudeBisectionRootFinder_->resetBoundaries( 75.0e3, 175.0e3 );

    // Set root-finder function as the area below the acceleration curve
    double estimatedLowerAltitudeBound = altitudeBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                    boost::bind( &lowerAltitudeBisectionFunction, _1, estimatedApoapsisRadius, planetaryRadius,
                                 planetaryGravitationalParameter, atmosphericInterfaceRadius, maximumAllowedHeatRate_,
                                 maximumAllowedHeatLoad_, densityFunction ) ) );
    std::cout << "Lower boundary: " << estimatedLowerAltitudeBound / 1.0e3 << " km" << std::endl;

    // Set root-finder boundaries for the upper altitude limits
    altitudeBisectionRootFinder_->resetBoundaries( 75.0e3, 175.0e3 );

    // Set root-finder function as the area below the acceleration curve
    double estimatedUpperAltitudeBound = altitudeBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                    boost::bind( &upperAltitudeBisectionFunction, _1, estimatedApoapsisRadius, planetaryRadius,
                                 planetaryGravitationalParameter, minimumAllowedDynamicPressure_, densityFunction ) ) );
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
    std::cout << "Tagert periapsis altitude: " << estimatedTargetPeriapsisAltitude / 1.0e3 << " km" << std::endl;

    // Save periapsis corridor altitudes to history
    historyOfEstimatedPeriapsisCorridorBoundaries_.push_back( std::make_pair( estimatedLowerAltitudeBound, estimatedUpperAltitudeBound ) );

    // Store corridor information to pair
    pairOfPeriapsisTargetingInformation_ = std::make_pair( isAltitudeWithinBounds, estimatedTargetPeriapsisAltitude );
}

} // namespace guidance_navigation_control

} // namespace tudat
