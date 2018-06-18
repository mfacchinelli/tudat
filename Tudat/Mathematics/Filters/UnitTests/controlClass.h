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

#ifndef TUDAT_CONTROL_CLASS_H
#define TUDAT_CONTROL_CLASS_H

#include <Eigen/Core>
#include <boost/function.hpp>

namespace tudat
{

namespace filters
{

//! Class for control vector
template< typename IndependentVariableType, typename DependentVariableType, int NumberOfElements >
class ControlWrapper
{
public:

    //! Typedef of the control vector.
    typedef Eigen::Matrix< DependentVariableType, NumberOfElements, 1 > DependentVector;

    //! Typedef of the function describing the system.
    typedef boost::function< DependentVector( const IndependentVariableType,
                                              const DependentVector& ) > ControlFunction;

    //! Default constructor.
    ControlWrapper( const ControlFunction& controlFunction ) :
        controlFunction_( controlFunction )
    { }

    //! Default destructor.
    ~ControlWrapper( ){ }

    DependentVector getControlVector( )
    {
        return controlVector_;
    }

    void setControlVector( const IndependentVariableType currentTime,
                           const DependentVector& currentStateVector )
    {
        controlVector_ = controlFunction_( currentTime, currentStateVector );
    }

private:

    ControlFunction controlFunction_;

    DependentVector controlVector_;

};


} // namespace filters

} // namespace tudat

#endif // TUDAT_CONTROL_CLASS_H