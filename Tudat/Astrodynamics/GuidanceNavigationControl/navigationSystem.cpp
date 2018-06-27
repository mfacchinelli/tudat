#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigationSystem.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Periapse time estimator (PTE).
void NavigationSystem::runPeriapseTimeEstimator( const std::map< double, Eigen::Vector3d >& estimatedAerodynamicAcceleration )
{
    // Extract aerodynamic accelerations of when the spacecraft is below the atmospheric interface altitude
    std::map< double, Eigen::Vector3d > estimatedAerodynamicAccelerationBelowAtmosphericInterface;
    for ( std::map< double, std::pair< Eigen::VectorXd, Eigen::VectorXd > >::const_iterator stateIterator =
          currentOrbitHistoryOfEstimatedStates_.begin( ); stateIterator != currentOrbitHistoryOfEstimatedStates_.end( );
          stateIterator++ )
    {
        if ( stateIterator->second.first.segment( 0, 3 ).norm( ) <= atmosphericInterfaceRadius_ )
        {
            estimatedAerodynamicAccelerationBelowAtmosphericInterface[ stateIterator->first ] =
                    estimatedAerodynamicAcceleration[ stateIterator->first ];
        }
    }


}

} // namespace navigation

} // namespace tudat
