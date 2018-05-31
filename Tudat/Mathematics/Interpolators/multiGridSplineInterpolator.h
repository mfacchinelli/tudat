/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *     Stackoverflow. C++ GCC4.4 warning: array subscript is above array bounds, 2009,
 *         http://stackoverflow.com/questions/
 *             1168525/c-gcc4-4-warning-array-subscript-is-above-array-bounds,
 *         last accessed: 27th December, 2013.
 *     GCC Mailing List. Re: How to fix 'array subscript is above array bounds' ?, 2012,
 *         http://gcc.gnu.org/ml/gcc-help/2012-04/msg00047.html, last accessed: 27th December,
 *         2013.
 *
 *    Notes
 *     Under older GCC-based compilers (4.3 and 4.4 series), it is known that this file will
 *     generate a spurious warning stating "array subscript is above array bounds". This warning
 *     can be safely ignored. It is recommended that you working with a GCC 4.5+ compiler, since
 *     the problem has been fixed in all versions that postdate 4.5. This warning has specifically
 *     been noted when compiling using the MinGW GCC 4.4.0 compiler under MS Windows. For more
 *     information on the nature of this warning, please take a look at Stackoverflow (2009) and
 *     GCC Mailing List (2012).
 *
 */

#ifndef TUDAT_MULTI_GRID_SPLINE_INTERPOLATOR_H
#define TUDAT_MULTI_GRID_SPLINE_INTERPOLATOR_H

#include "Tudat/Mathematics/BasicMathematics/nearestNeighbourSearch.h"
#include "Tudat/Mathematics/Interpolators/multiDimensionalInterpolator.h"

namespace tudat
{

namespace interpolators
{

//! Class for performing multi-dimensional grid spline interpolation for an arbitrary number of
//! independent variables.
/*!
 *  Class for performing multi-dimensional grid spline interpolation for an arbitrary number of
 *  independent variables. Note that the types (i.e. double, float) of all independent variables
 *  must be the same.
 *  \tparam IndependentVariableType Type for independent variables.
 *  \tparam DependentVariableType Type for dependent variable.
 *  \tparam NumberOfDimensions Number of independent variables.
 */
template< typename IndependentVariableType, typename DependentVariableType, int NumberOfDimensions >
class MultiGridSplineInterpolator: public MultiDimensionalInterpolator< IndependentVariableType,
        DependentVariableType, NumberOfDimensions >
{
public:

    // Using statements to prevent having to put 'this' everywhere in the code.
    using MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >::
    dependentData_;
    using MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >::
    independentValues_;
    using MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >::
    lookUpSchemes_;
    using MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >::
    defaultExtrapolationValue_;

    //! Default constructor taking independent and dependent variable data.
    /*!
     *  Default constructor taking independent and dependent variable data.
     *  \param independentValues Vector of vectors containing data points of independent variables,
     *      each must be sorted in ascending order.
     *  \param dependentData Multi-dimensional array of dependent data at each point of
     *      hyper-rectangular grid formed by independent variable points.
     *  \param selectedLookupScheme Identifier of lookupscheme from enum. This algorithm is used
     *      to find the nearest lower data point in the independent variables when requesting
     *      interpolation.
     *  \param boundaryHandling Vector of boundary handling methods, in case independent variable is outside the
     *      specified range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    MultiGridSplineInterpolator( const std::vector< std::vector< IndependentVariableType > > independentValues,
                                 const boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) >
                                 dependentData,
                                 const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                                 const std::vector< BoundaryInterpolationType > boundaryHandling =
            std::vector< BoundaryInterpolationType >( NumberOfDimensions, extrapolate_at_boundary ),
                                 const DependentVariableType defaultExtrapolationValue =
            IdentityElement< DependentVariableType >::getAdditionIdentity( ) ) :
        MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions >(
            boundaryHandling, defaultExtrapolationValue )
    {
        // Save (in)dependent variables
        independentValues_ = independentValues;
        dependentData_.resize( reinterpret_cast< boost::array< size_t,
                               boost::multi_array< DependentVariableType,
                               static_cast< size_t >( NumberOfDimensions ) >::dimensionality > const& >(
                                   *dependentData.shape( ) ) ); // resize dependent data container
        dependentData_ = dependentData;

        // Check consistency of template arguments and input variables.
        if ( independentValues.size( ) != NumberOfDimensions )
        {
            throw std::runtime_error( "Error: dimension of independent value vector provided to constructor "
                                      "incompatible with template parameter." );
        }

        // Check consistency of input data of dependent and independent data.
        for ( int i = 0; i < NumberOfDimensions; i++ )
        {
            if ( independentValues[ i ].size( ) != dependentData.shape( )[ i ] )
            {
                std::string errorMessage = "Error: number of data points in dimension " +
                        std::to_string( i ) + " of independent and dependent data incompatible.";
                throw std::runtime_error( errorMessage );
            }
        }

        // Create lookup scheme from independent variable data points.
        this->makeLookupSchemes( selectedLookupScheme );
    }

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~MultiGridSplineInterpolator( ){ }

    //! Function to perform interpolation.
    /*!
     *  This function performs the multilinear interpolation.
     *  \param independentValuesToInterpolate Vector of values of independent variables at which
     *      the value of the dependent variable is to be determined.
     *  \return Interpolated value of dependent variable in all dimensions.
     */
    DependentVariableType interpolate( const std::vector< IndependentVariableType >& independentValuesToInterpolate )
    {
        // Check whether size of independent variable vector is correct
        if ( independentValuesToInterpolate.size( ) != NumberOfDimensions )
        {
            throw std::runtime_error( "Error in multi-dimensional interpolator. The number of independent variables "
                                      "provided is incompatible with the previous definition. Provided: " +
                                      std::to_string( independentValuesToInterpolate.size( ) ) + ". Needed: " +
                                      std::to_string( NumberOfDimensions ) );
        }

        // Create copy of values to interpolate, such that it can be modified
        std::vector< IndependentVariableType > localIndependentValuesToInterpolate = independentValuesToInterpolate;

        // Check that independent variables are in range
        bool useDefault = false;
        for ( unsigned int i = 0; i < NumberOfDimensions; i++ )
        {
            this->checkBoundaryCase( localIndependentValuesToInterpolate.at( i ), i, useDefault );
            if ( useDefault )
            {
                return defaultExtrapolationValue_;
            }
        }

    }

private:



};

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_MULTI_GRID_SPLINE_INTERPOLATOR_H
