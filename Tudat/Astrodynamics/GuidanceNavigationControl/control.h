/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef MICHELE_GNC_CONTROL
#define MICHELE_GNC_CONTROL

namespace tudat
{

namespace guidance_navigation_control
{

//! Class for control system of an aerobraking maneuver.
class ControlSystem
{
public:

    //! Constructor.
    ControlSystem( ) { }

    //! Destructor.
    ~ControlSystem( ) { }

    Eigen::VectorXd getControlVector( )
    {
        return controlVector_;
    }

private:

    Eigen::VectorXd controlVector_;

};

} // namespace thesis

} // namespace tudat

#endif // MICHELE_GNC_CONTROL
