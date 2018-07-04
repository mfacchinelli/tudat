/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_CONTROL_TORQUE_H
#define TUDAT_CONTROL_TORQUE_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/GuidanceNavigationControl/controlSystem.h"

namespace tudat
{

namespace guidance_navigation_control
{

//! Class to link the control system to the torque model.
/*!
 *
 */
class ControlTorque : public basic_astrodynamics::TorqueModel
{
public:

    //! Constructor
    ControlTorque( const boost::shared_ptr< ControlSystem > controlSystem ) :
        controlSystem_( controlSystem )
    { }

    //! Destructor
    ~ControlTorque( ) { }

    //! Function to retrieve the current value of the torque.
    /*!
     *  Function to retrieve the current value of the torque, which is computed by the
     *  control system based on the current estimated state.
     *  \return Current torque, as computed by last call to updateMembers function.
     */
    Eigen::Vector3d getTorque( )
    {
        // Retrieve and return the control vector from the control system.
        return controlSystem_->getCurrentAttitudeControlVector( );
    }

    //! Update member variables used by the torque model.
    /*!
     *  Update member variables used by the torque model. This function is empty, since the
     *  control system is automatically updated by the onboard computer object.
     */
    void updateMembers( const double currentTime )
    {
        // nothing to be done, the control system is updated by the onboard computer
    }

protected:

private:

    //! Pointer to the control system to be used to retrieve the torque.
    const boost::shared_ptr< ControlSystem > controlSystem_;

};

} // namespace guidance_navigation_control

} // namespace tudat

#endif // TUDAT_CONTROL_TORQUE_H
