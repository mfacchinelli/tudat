#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigation.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Periapse time estimator (PTE).
void NavigationSystem::runPeriapseTimeEstimator( std::map< double, std::pair< Eigen::VectorXd, Eigen::VectorXd > >& estimatedState,
                                              const std::map< double, Eigen::Vector3d >& estimatedAerodynamicAcceleration )
{

}

} // namespace navigation

} // namespace tudat
