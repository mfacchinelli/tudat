#include "Tudat/Astrodynamics/GuidanceNavigationControl/navigation.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Periapse time estimator (PTE).
void NavigationSystem::periapseTimeEstimator( std::map< double, std::pair< Eigen::VectorXd, Eigen::VectorXd > >& estimatedState,
                                              const std::map< double, Eigen::Vector3d >& estimatedAerodynamicAcceleration,
                                              const boost::function< Eigen::MatrixXd( const double, const Eigen::VectorXd& ) >& stateTransitionMatrixFunction )
{

}

} // namespace navigation

} // namespace tudat
