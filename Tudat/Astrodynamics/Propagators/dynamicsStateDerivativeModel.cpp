#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/thirdBodyGravityPartial.h"

namespace tudat
{
namespace propagators
{

//! Function to retrieve specific acceleration partial object from list of state derivative partials
boost::shared_ptr< acceleration_partials::AccelerationPartial > getAccelerationPartialForBody(
        const orbit_determination::StateDerivativePartialsMap& accelerationPartials,
        const basic_astrodynamics::AvailableAcceleration accelerationType,
        const std::string& bodyExertingAcceleration,
        const std::string& bodyUndergoignAcceleration,
        const std::string& centralBody )
{
    boost::shared_ptr< acceleration_partials::AccelerationPartial > matchingAccelerationPartial;

    // Iterate over all partials
    for( unsigned int i = 0; i < accelerationPartials.size( ); i++ )
    {
        for( unsigned int j = 0; j < accelerationPartials.at( i ).size( ); j++ )
        {
            // Check id state derivative type is acceleration
            boost::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial =
                    boost::dynamic_pointer_cast< acceleration_partials::AccelerationPartial >(
                        accelerationPartials.at( i ).at( j ) );
            if( accelerationPartial != NULL )
            {
                // Compare current partial against required data
                bool partialIdentified = false;
                if( ( accelerationPartial->getAccelerationType( ) == accelerationType ) &&
                        ( accelerationPartial->getAcceleratedBody( ) == bodyExertingAcceleration ) &&
                        ( accelerationPartial->getAcceleratingBody( ) == bodyUndergoignAcceleration ) )
                {
                    if( !basic_astrodynamics::isAccelerationFromThirdBody( accelerationType  ) )
                    {
                        partialIdentified = true;
                    }
                    else if( getCentralBodyNameFromThirdBodyAccelerationPartial( accelerationPartial ) == centralBody )
                    {
                        partialIdentified = true;
                    }
                }

                // Set output; check if multiple models are identified.
                if( partialIdentified )
                {
                    if( matchingAccelerationPartial != NULL )
                    {
                        throw std::runtime_error( "Error when getting acceleration partial, found multiple matching accelerations." );
                    }
                    else
                    {
                        matchingAccelerationPartial = accelerationPartial;
                    }
                }
            }
        }
    }

    return matchingAccelerationPartial;
}

} // namespace propagators

} // namespace tudat
