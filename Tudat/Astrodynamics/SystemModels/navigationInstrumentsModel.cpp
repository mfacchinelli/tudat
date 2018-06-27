#include "Tudat/Astrodynamics/SystemModels/navigationInstrumentsModel.h"

namespace tudat
{

namespace system_models
{

//! Function to compute the scale-misalignment matrix from the scale and misalignment vectors.
Eigen::Matrix3d computeScaleMisalignmentMatrix( const Eigen::Vector3d& scaleFactorVector,
                                                const Eigen::Vector6d& misalignmentVector )
{
    // Compute misalignment matrix
    Eigen::Matrix3d misalignmentMatrix = Eigen::Matrix3d::Zero( );
    misalignmentMatrix( 0, 1 ) = misalignmentVector[ 0 ];
    misalignmentMatrix( 0, 2 ) = misalignmentVector[ 1 ];
    misalignmentMatrix( 1, 0 ) = misalignmentVector[ 2 ];
    misalignmentMatrix( 1, 2 ) = misalignmentVector[ 3 ];
    misalignmentMatrix( 2, 0 ) = misalignmentVector[ 4 ];
    misalignmentMatrix( 2, 1 ) = misalignmentVector[ 5 ];

    // Give output
    return Eigen::Matrix3d::Identity( ) + scaleFactorVector.asDiagonal( ) + misalignmentMatrix;
}

} // namespace system_models

} // namespace tudat
