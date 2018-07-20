#include "Tudat/Astrodynamics/GuidanceNavigationControl/guidanceSystem.h"

#include "Tudat/Mathematics/BasicMathematics/functionProxy.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/extraFunctions.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Function to run corridor estimator (CE).
void GuidanceSystem::runCorridorEstimator( const unsigned int currentOrbitCounter )
{
    // Compute predicted apoapsis altitude
    double predictedApoapsisAltitude;

    // Compute predicted periapsis altitude
    double predictedPeriapsisAltitude;

    // Set root-finder boundaries for the lower altitude limits
    altitudeBisectionRootFinder_->resetBoundaries( 75.0e3, 125.0e3 );

    // Set root-finder function as the area below the acceleration curve
    double estimatedLowerAltitudeBound = altitudeBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                    boost::bind( &lowerAltitudeBisectionFunction, _1,  ) ) );

    // Set root-finder boundaries for the upper altitude limits
    altitudeBisectionRootFinder_->resetBoundaries( 100.0e3, 150.0e3 );

    // Set root-finder function as the area below the acceleration curve
    double estimatedUpperAltitudeBound = altitudeBisectionRootFinder_->execute(
                boost::make_shared< basic_mathematics::FunctionProxy< double, double > >(
                    boost::bind( &upperAltitudeBisectionFunction, _1,  ) ) );

    // Check whether predicted altitude is within bounds
    bool isAltitudeWithinBounds = ( ( predictedPeriapsisAltitude > estimatedLowerAltitudeBound ) &&
                                    ( predictedPeriapsisAltitude < estimatedUpperAltitudeBound ) );

    // Compute target altitude
    // If the predicted periapsis altitude falls within the bounds, then the target altitude is simply the
    // predicted one, otherwise it is the mid-point of the corridor
    double estimatedTargetPeriapsisAltitude = isAltitudeWithinBounds ? predictedPeriapsisAltitude :
                                                                       0.5 * ( estimatedLowerAltitudeBound +
                                                                               estimatedUpperAltitudeBound );

    // Save periapsis corridor altitudes to history
    historyOfEstimatedPeriapsisCorridorBoundaries_[ currentOrbitCounter ] = std::make_pair( estimatedLowerAltitudeBound,
                                                                                            estimatedUpperAltitudeBound );

    // Store corridor information to pair
    pairOfPeriapsisTargetingInformation_ = std::make_pair( isAltitudeWithinBounds, estimatedTargetPeriapsisAltitude );
}

} // namespace guidance_navigation_control

} // namespace tudat
