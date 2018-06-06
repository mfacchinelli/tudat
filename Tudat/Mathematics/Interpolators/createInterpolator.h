/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEINTERPOLATOR_H
#define TUDAT_CREATEINTERPOLATOR_H

#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/hermiteCubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Mathematics/Interpolators/piecewiseConstantInterpolator.h"

#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"

#include "Tudat/InputOutput/mapTextFileReader.h"

namespace tudat
{

namespace interpolators
{

//! Enum of available interpolator types.
enum InterpolatorTypes
{
    linear_interpolator = 0,
    multi_linear_interpolator = 1,
    cubic_spline_interpolator = 2,
    lagrange_interpolator = 3,
    hermite_spline_interpolator = 4,
    piecewise_constant_interpolator = 5
};

//! Base class for providing settings for creating an interpolator.
/*!
 *  Base class for providing settings for creating an interpolator using the createInterpolator
 *  function. This base class is not-functional, i.e. a derived class needs to be used.
 */
template< unsigned int NumberOfDimensions >
class InterpolatorSettings
{  
public:

    //! Default constructor.
    /*!
     *  Default constructor. Constructor taking a vector of boundary handling methods. The vector length needs
     *  to be equal to the number of dimensions.
     *  \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     *  \param boundaryHandling Vector of boundary handling methods, in case independent variable is outside the
     *      specified range.
     */
    InterpolatorSettings( const InterpolatorTypes interpolatorType,
                          const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                          const std::vector< BoundaryInterpolationType >& boundaryHandling =
            std::vector< BoundaryInterpolationType >( NumberOfDimensions, extrapolate_at_boundary ) ) :
        interpolatorType_( interpolatorType ), selectedLookupScheme_( selectedLookupScheme ),
        boundaryHandling_( boundaryHandling )
    {
        // Check if dimensions are consistent
        if ( NumberOfDimensions != boundaryHandling.size( ) )
        {
            throw std::runtime_error( "Error while creating interpolator settings. Inconsistent dimensions in "
                                      "boundary handling vector." );
        }

        // Check that if interpolator type matches with number of dimensions
        std::vector< bool > isMethodOneDimensional = std::vector< bool >( 6, true );
        isMethodOneDimensional.at( static_cast< unsigned int >( multi_linear_interpolator ) ) = false;
        if ( NumberOfDimensions > 1 && isMethodOneDimensional.at( static_cast< unsigned int >( interpolatorType_ ) ) )
        {
            throw std::runtime_error( "Error while creating interpolator settings. Number of dimensions is greater "
                                      "than 1, but a one-dimensional interpolator has been selected." );
        }
    }

    //! Constructor.
    /*!
     *  Constructor.
     *  \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     *  \param boundaryHandling Boundary handling method, in case independent variable is outside the
     *      specified range.
     */
    InterpolatorSettings( const InterpolatorTypes interpolatorType,
                          const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
                          const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary ) :
        InterpolatorSettings( interpolatorType, selectedLookupScheme,
                              std::vector< BoundaryInterpolationType >( NumberOfDimensions, boundaryHandling ) )
    { }

    //! Virtual destructor.
    virtual ~InterpolatorSettings( ) { }

    //! Function to get the selected type of interpolator.
    /*!
     * Function to get the selected type of interpolator.
     * \return Selected type of interpolator.
     */
    InterpolatorTypes getInterpolatorType( )
    {
        return interpolatorType_;
    }

    //! Function to get the selected type of lookup scheme for independent variables.
    /*!
     * Function to get the selected type of lookup scheme for independent variables.
     * \return Selected type of lookup scheme for independent variables.
     */
    AvailableLookupScheme getSelectedLookupScheme( )
    {
        return selectedLookupScheme_;
    }

    //! Function to retrieve boundary handling method.
    /*!
     * Function to retrieve boundary handling method.
     * \return Boundary handling method.
     */
    std::vector< BoundaryInterpolationType > getBoundaryInterpolationType( )
    {
        return boundaryHandling_;
    }

protected:

    //! Selected type of interpolator.
    InterpolatorTypes interpolatorType_;

    //! Selected type of lookup scheme for independent variables.
    AvailableLookupScheme selectedLookupScheme_;

    //! Boundary handling method.
    std::vector< BoundaryInterpolationType > boundaryHandling_;

};

//! Base class for providing settings for creating a one-dimensional interpolator.
/*!
 *  Base class for providing settings for creating an interpolator using the createInterpolator
 *  function. This base class is functional, i.e. may be used directly for interpolator types
 *  requiring no additional information. For interpolators that do require more information,
 *  a derived class is provided.
 */
class OneDimensionalInterpolatorSettings : public InterpolatorSettings< 1 >
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     * \param interpolatorType Selected type of interpolator.
     * \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     * \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double,
     *      time step is a double if false.
     * \param boundaryHandling Boundary handling method in case independent variable is outside the
     *      specified range.
     * \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    OneDimensionalInterpolatorSettings(
            const InterpolatorTypes interpolatorType,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
            const bool useLongDoubleTimeStep = false,
            const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary ) :
        InterpolatorSettings< 1 >( interpolatorType, selectedLookupScheme, boundaryHandling ),
        useLongDoubleTimeStep_( useLongDoubleTimeStep )
    { }

    //! Virtual destructor
    virtual ~OneDimensionalInterpolatorSettings( ){ }

    //! Function to reset the use of long double type for time step.
    /*!
     * Function to reset the use of long double type for time step.
     * \param Boolean denoting whether time step is to be a long double.
     */
    void resetUseLongDoubleTimeStep( const bool useLongDoubleTimeStep )
    {
        useLongDoubleTimeStep_ = useLongDoubleTimeStep;
    }

    //! Function to get a boolean denoting whether time step is to be a long double.
    /*!
     * Function to get a boolean denoting whether time step is to be a long double.
     * \return Boolean denoting whether time step is to be a long double.
     */
    bool getUseLongDoubleTimeStep( )
    {
        return useLongDoubleTimeStep_;
    }

protected:

    //! Boolean denoting whether time step is to be a long double.
    bool useLongDoubleTimeStep_;

};

//! Class for providing settings to creating a Lagrange interpolator.
class LagrangeInterpolatorSettings : public OneDimensionalInterpolatorSettings
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param interpolatorOrder Order of the Lagrange interpolator that is to be created.
     * \param useLongDoubleTimeStep Boolean denoting whether time step is to be a long double,
     * time step is a double if false.
     * \param selectedLookupScheme Selected type of lookup scheme for independent variables.
     * \param boundaryHandling Variable denoting the method by which the boundary interpolation is handled.
     */
    LagrangeInterpolatorSettings(
            const int interpolatorOrder,
            const bool useLongDoubleTimeStep = 0,
            const AvailableLookupScheme selectedLookupScheme = huntingAlgorithm,
            const LagrangeInterpolatorBoundaryHandling lagrangeBoundaryHandling = lagrange_cubic_spline_boundary_interpolation,
            const BoundaryInterpolationType boundaryHandling = extrapolate_at_boundary ) :
        OneDimensionalInterpolatorSettings( lagrange_interpolator, selectedLookupScheme, useLongDoubleTimeStep, boundaryHandling ),
        interpolatorOrder_( interpolatorOrder ),
        lagrangeBoundaryHandling_( lagrangeBoundaryHandling )
    { }

    //! Destructor
    ~LagrangeInterpolatorSettings( ){ }

    //! Function to get the order of the Lagrange interpolator that is to be created.
    /*!
     * Function to get the order of the Lagrange interpolator that is to be created.
     * \return Order of the Lagrange interpolator that is to be created.
     */
    int getInterpolatorOrder( )
    {
        return interpolatorOrder_;
    }

    LagrangeInterpolatorBoundaryHandling getBoundaryHandling( )
    {
        return lagrangeBoundaryHandling_;
    }


protected:

    //! Order of the Lagrange interpolator that is to be created.
    int interpolatorOrder_;

    LagrangeInterpolatorBoundaryHandling lagrangeBoundaryHandling_;

};


//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief DataMapSettings
 * The class can be used to provided a map directly, or as a base class for loading the map using different procedures.
 */
template< typename IndependentType, typename DependentType >
class DataMapSettings
{
public:

    //! Empty constructor.
    /*!
     * Empty constructor.
     */
    DataMapSettings( ) { }

    //! Constructor with a data map.
    /*!
     * Constructor to be used when the data map is provided directly.
     * \param dataMap The data map containing values for the independent and dependent variables.
     */
    DataMapSettings( const std::map< IndependentType, DependentType >& dataMap ) : dataMap_( dataMap ) { }

    //! Virtual destructor
    virtual ~DataMapSettings( ) { }

    //! Get the data map associated to the current setting object.
    /*!
     * @copybrief getDataMap
     * \return The data map associated to the current setting object.
     */
    virtual std::map< IndependentType, DependentType > getDataMap( ) const
    {
        return dataMap_;
    }

protected:

    //! The data map directly provided by the user in the constructor.
    const std::map< IndependentType, DependentType > dataMap_;

};

//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief IndependentDependentDataMapSettings
 * The data map will be created by combining values of independent and dependent variables.
 */
template< typename IndependentType, typename DependentType >
class IndependentDependentDataMapSettings : public DataMapSettings< IndependentType, DependentType >
{
public:

    //! Consturctor.
    /*!
     * Constructor.
     * \param independentVariableValues Vector containing the values of the indepedent variable.
     * \param dependentVariableValues Vector containing the values of the depedent variable.
     */
    IndependentDependentDataMapSettings( const std::vector< IndependentType >& independentVariableValues,
                                         const std::vector< DependentType >& dependentVariableValues ) :
        DataMapSettings< IndependentType, DependentType >( ),
        independentVariableValues_( independentVariableValues ),
        dependentVariableValues_( dependentVariableValues ) { }

    //! Virtual destructor
    virtual ~IndependentDependentDataMapSettings( ) { }

    //! Vector containing the values of the indepedent variable.
    std::vector< IndependentType > independentVariableValues_;

    //! Vector containing the values of the depedent variable.
    std::vector< DependentType > dependentVariableValues_;

    //! Get the data map associated to the current setting object.
    /*!
     * @copybrief getDataMap
     * The map is created from the provided vectors of independent and depedent variables values.
     * \return The data map associated to the current setting object.
     */
    virtual std::map< IndependentType, DependentType > getDataMap( ) const
    {
        if ( independentVariableValues_.size( ) != dependentVariableValues_.size( ) )
        {
            std::cerr << "Could not get data map because the size of the independent and dependent variables values "
                         "is inconsistent." << std::endl;
            throw;
        }
        std::map< IndependentType, DependentType > dataMap;
        for ( unsigned int i = 0; i < independentVariableValues_.size( ); ++i )
        {
            dataMap[ independentVariableValues_.at( i ) ] = dependentVariableValues_.at( i );
        }
        return dataMap;
    }

};

//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief FromFileDataMapSettings
 * The data map will be read from a file.
 */
template< typename EigenVectorType >
class FromFileDataMapSettings : public DataMapSettings< typename EigenVectorType::Scalar, EigenVectorType >
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param relativeFilePath Relative path to the file from which the map is to be loaded.
     */
    FromFileDataMapSettings( const std::string& relativeFilePath ) :
        DataMapSettings< typename EigenVectorType::Scalar, EigenVectorType >( ),
        relativeFilePath_( relativeFilePath ) { }

    //! Virtual destructor
    virtual ~FromFileDataMapSettings( ) { }

    //! Relative path to the file from which the map is to be loaded.
    std::string relativeFilePath_;

    //! Get the data map associated to the current setting object.
    /*!
     * @copybrief getDataMap
     * The map is loaded from the specified relativeFilePath_.
     * \return The data map associated to the current setting object.
     */
    virtual std::map< typename EigenVectorType::Scalar, EigenVectorType > getDataMap( ) const
    {
        return input_output::readEigenVectorMapFromFile< EigenVectorType >( relativeFilePath_ );
    }

};

//! Class defening the settings to be used to create a map of data (used for interpolation).
/*!
 * @copybrief HermiteDataSettings
 * The data map will be created directly, in addition to the first derivatives.
 */
template< typename IndependentType, typename DependentType >
class HermiteDataSettings : public DataMapSettings< IndependentType, DependentType >
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param dataToInterpolate The data map containing values for the independent and dependent variables.
     * \param firstDerivativeOfDependentVariables Vector containing the first derivatives of the depedent variables.
     */
    HermiteDataSettings( const std::map< IndependentType, DependentType >& dataToInterpolate,
                         const std::vector< DependentType >& firstDerivativeOfDependentVariables ) :
        DataMapSettings< IndependentType, DependentType >( dataToInterpolate ),
        firstDerivativeOfDependentVariables_( firstDerivativeOfDependentVariables ) { }

    //! Virtual destructor
    virtual ~HermiteDataSettings( ) { }

    //! Vector containing the first derivatives of the depedent variables.
    std::vector< DependentType > firstDerivativeOfDependentVariables_;

};

//! Class containing (the settings to create) the data needed for the interpolation and the settings to create the
//! interpolator.
/*!
 * @copybrief DataInterpolationSettings
 */
template< typename IndependentType, typename DependentType >
class DataInterpolationSettings
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param dataSettings Object containing (the settings to create) the data needed for the interpolation.
     * \param interpolatorSettings Object containing the settings to create the interpolator to be used.
     */
    DataInterpolationSettings(
            const boost::shared_ptr< DataMapSettings< IndependentType, DependentType > >& dataSettings,
            const boost::shared_ptr< OneDimensionalInterpolatorSettings >& interpolatorSettings ) :
        dataSettings_( dataSettings ), interpolatorSettings_( interpolatorSettings ) { }

    //! Virtual destructor
    virtual ~DataInterpolationSettings( ){ }

    //! Object containing (the settings to create) the data needed for the interpolation.
    boost::shared_ptr< DataMapSettings< IndependentType, DependentType > > dataSettings_;

    //! Object containing the settings to create the interpolator to be used.
    boost::shared_ptr< OneDimensionalInterpolatorSettings > interpolatorSettings_;

};

//! Function to create a one-dimensional interpolator
/*!
 *  Function to create a one-dimensional interpolator from the data that is to be interpolated,
 *  as well as the settings that are to be used to create the interpolator.
 *  \param dataToInterpolate Map providing data that is to be interpolated (key = independent
 *  variables, value = dependent variables)
 *  \param interpolatorSettings Settings that are to be used to create interpolator
 *  \param firstDerivativeOfDependentVariables First derivative of dependent variables w.r.t. independent variable at
 *  independent variables values in values of dataToInterpolate. By default, this vector is empty, it only needs to
 *  be supplied if the selected interpolator requires this data (e.g. Hermite spline).
 *  \return Interpolator created from dataToInterpolate using interpolatorSettings.
 */
template< typename IndependentVariableType, typename DependentVariableType >
boost::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
createOneDimensionalInterpolator(
        const std::map< IndependentVariableType, DependentVariableType > dataToInterpolate,
        const boost::shared_ptr< OneDimensionalInterpolatorSettings > interpolatorSettings,
        const std::vector< DependentVariableType > firstDerivativeOfDependentVariables =
        std::vector< DependentVariableType >( ) )
{
    boost::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
            createdInterpolator;

    // Check type of interpolator.
    switch( interpolatorSettings->getInterpolatorType( ) )
    {
    case linear_interpolator:
        createdInterpolator = boost::make_shared< LinearInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryInterpolationType( ).at( 0 ) );
        break;
    case cubic_spline_interpolator:
    {
        if( !interpolatorSettings->getUseLongDoubleTimeStep( ) )
        {
            createdInterpolator = boost::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType > >(
                        dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                        interpolatorSettings->getBoundaryInterpolationType( ).at( 0 ) );
        }
        else
        {
            createdInterpolator = boost::make_shared< CubicSplineInterpolator
                    < IndependentVariableType, DependentVariableType, long double > >(
                        dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                        interpolatorSettings->getBoundaryInterpolationType( ).at( 0 ) );
        }
        break;
    }
    case lagrange_interpolator:
    {
        // Check consistency of input
        boost::shared_ptr< LagrangeInterpolatorSettings > lagrangeInterpolatorSettings =
                boost::dynamic_pointer_cast< LagrangeInterpolatorSettings >( interpolatorSettings );
        if( lagrangeInterpolatorSettings != NULL )
        {
            // Create Lagrange interpolator with requested time step type
            if( !lagrangeInterpolatorSettings->getUseLongDoubleTimeStep( ) )
            {
                createdInterpolator = boost::make_shared< LagrangeInterpolator
                        < IndependentVariableType, DependentVariableType, double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ),
                            lagrangeInterpolatorSettings->getBoundaryHandling( ),
                            interpolatorSettings->getBoundaryInterpolationType( ).at( 0 ) );
            }
            else
            {
                createdInterpolator = boost::make_shared< LagrangeInterpolator
                        < IndependentVariableType, DependentVariableType, long double > >(
                            dataToInterpolate, lagrangeInterpolatorSettings->getInterpolatorOrder( ),
                            interpolatorSettings->getSelectedLookupScheme( ),
                            lagrangeInterpolatorSettings->getBoundaryHandling( ),
                            interpolatorSettings->getBoundaryInterpolationType( ).at( 0 ) );
            }
        }
        else
        {
            throw std::runtime_error( "Error, did not recognize lagrange interpolator settings" );

        }
        break;
    }
    case hermite_spline_interpolator:
    {
        if( firstDerivativeOfDependentVariables.size( ) != dataToInterpolate.size( ) )
        {
            throw std::runtime_error(
                        "Error when creating hermite spline interpolator, derivative size is inconsistent" );
        }
        createdInterpolator = boost::make_shared< HermiteCubicSplineInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, firstDerivativeOfDependentVariables,
                    interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryInterpolationType( ).at( 0 ) );
        break;
    }
    case piecewise_constant_interpolator:
        createdInterpolator = boost::make_shared< PiecewiseConstantInterpolator
                < IndependentVariableType, DependentVariableType > >(
                    dataToInterpolate, interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryInterpolationType( ).at( 0 ) );
        break;
    default:
        throw std::runtime_error( "Error when making interpolator, function cannot be used to create interplator of type " +
                                  std::to_string( interpolatorSettings->getInterpolatorType( ) ) );
    }
    return createdInterpolator;
}

//! Function to create an interpolator from DataInterpolationSettings
/*!
 *  Function to create an interpolator from DataInterpolationSettings
 *  \param dataInterpolationSettings Object containing the data that is to be interpolated and settings that are to be
 *  used to create the interpolator.
 *  \return Interpolator created from dataToInterpolate using interpolatorSettings.
 */
template< typename IndependentType, typename DependentType >
boost::shared_ptr< OneDimensionalInterpolator< IndependentType, DependentType > > createOneDimensionalInterpolator(
        const boost::shared_ptr< DataInterpolationSettings< IndependentType, DependentType > >
        dataInterpolationSettings )
{
    std::vector< DependentType > firstDerivativeOfDependentVariables;
    boost::shared_ptr< HermiteDataSettings< IndependentType, DependentType > > hermiteDataSettings =
            boost::dynamic_pointer_cast< HermiteDataSettings< IndependentType, DependentType > >(
                dataInterpolationSettings->dataSettings_ );
    if ( hermiteDataSettings )
    {
        firstDerivativeOfDependentVariables = hermiteDataSettings->firstDerivativeOfDependentVariables_;
    }
    return createOneDimensionalInterpolator( dataInterpolationSettings->dataSettings_->getDataMap( ),
                                             dataInterpolationSettings->interpolatorSettings_,
                                             firstDerivativeOfDependentVariables );
}

//! Function to create a multi-dimensional interpolator
/*!
 *  Function to create a multi-dimensional interpolator from the data that is to be interpolated,
 *  as well as the settings that are to be used to create the interpolator.
 *  \param independentValues Vector of vectors containing data points of independent variables,
 *  each must be sorted in ascending order.
 *  \param dependentData Multi-dimensional array of dependent data at each point of
 *  hyper-rectangular grid formed by independent variable points.
 *  \param interpolatorSettings Settings that are to be used to create interpolator
 *  \return Interpolator created from independentValues and dependentData using interpolatorSettings.
 */
template< typename IndependentVariableType, typename DependentVariableType, unsigned int NumberOfDimensions >
boost::shared_ptr< MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions > >
createMultiDimensionalInterpolator(
        const std::vector< std::vector< IndependentVariableType > >& independentValues,
        const boost::multi_array< DependentVariableType, static_cast< size_t >( NumberOfDimensions ) >& dependentData,
        const boost::shared_ptr< InterpolatorSettings< NumberOfDimensions > > interpolatorSettings )
{
    boost::shared_ptr< MultiDimensionalInterpolator< IndependentVariableType, DependentVariableType, NumberOfDimensions > >
            createdInterpolator;

    // Check type of interpolator.
    switch ( interpolatorSettings->getInterpolatorType( ) )
    {
    case multi_linear_interpolator:
        createdInterpolator = boost::make_shared< MultiLinearInterpolator
                < IndependentVariableType, DependentVariableType, NumberOfDimensions > >(
                    independentValues, dependentData, interpolatorSettings->getSelectedLookupScheme( ),
                    interpolatorSettings->getBoundaryInterpolationType( ) );
        break;
    default:
        throw std::runtime_error( "Error when making interpolator, function cannot be used to create interplator of type " +
                                  std::to_string( interpolatorSettings->getInterpolatorType( ) ) );
    }
    return createdInterpolator;
}

} // namespace interpolators

} // namespace tudat

#endif // TUDAT_CREATEINTERPOLATOR_H
