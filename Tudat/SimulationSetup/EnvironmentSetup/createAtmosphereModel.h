/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEATMOSPHEREMODEL_H
#define TUDAT_CREATEATMOSPHEREMODEL_H

#include <string>
#include <map>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/atmosphereModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Basics/identityElements.h"

namespace tudat
{

namespace simulation_setup
{

using namespace aerodynamics;

//! List of wind models available in simulations
/*!
 *  List of wind models available in simulations. Wind models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum WindModelTypes
{
    custom_wind_model
};

//! Class for providing settings for wind model.
/*!
 *  Class for providing settings for automatic wind model creation. This class is a
 *  functional (base) class for settings of wind models that require no information in
 *  addition to their type. Wind model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
class WindModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param windModelType Type of wind model that is to be created
     */
    WindModelSettings( const WindModelTypes windModelType ):
        windModelType_( windModelType ){ }

    //! Destructor
    virtual ~WindModelSettings( ){ }

    //! Function to retrieve type of wind model that is to be created
    /*!
     * Function to retrieve type of wind model that is to be created
     * \return Type of wind model that is to be created
     */
    WindModelTypes getWindModelType( )
    {
        return windModelType_;
    }

protected:

    //! Type of wind model that is to be created
    WindModelTypes windModelType_;
};

//! Class to define settings for a custom, user-defined, wind model
class CustomWindModelSettings: public WindModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     */
    CustomWindModelSettings(
            const boost::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction ):
        WindModelSettings( custom_wind_model ), windFunction_( windFunction ){ }

    //! Destructor
    ~CustomWindModelSettings( ){ }

    //! Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
    /*!
     * Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
     * \return Function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    boost::function< Eigen::Vector3d( const double, const double, const double, const double ) > getWindFunction( )
    {
        return windFunction_;
    }

    //! Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
    /*!
     * Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
     * \param windFunction New function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    void setWindFunction(
            const boost::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction )
    {
        windFunction_ = windFunction;
    }

protected:

    //! Function that returns wind vector as a function of altitude, longitude, latitude and time (in that order).
    boost::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;
};

//! List of atmosphere models available in simulations
/*!
 *  List of atmosphere models available in simulations. Atmosphere models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum AtmosphereTypes
{
    exponential_atmosphere,
    tabulated_atmosphere,
    nrlmsise00
};

//! Class for providing settings for atmosphere model.
/*!
 *  Class for providing settings for automatic atmosphere model creation. This class is a
 *  functional (base) class for settings of atmosphere models that require no information in
 *  addition to their type. Atmosphere model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
class AtmosphereSettings
{
public:

    //! Constructor, sets type of atmosphere model.
    /*!
     *  Constructor, sets type of atmosphere model. Settings for atmosphere models requiring
     *  additional information should be defined in a derived class.
     *  \param atmosphereType Type of atmosphere model that is to be created.
     */
    AtmosphereSettings( const AtmosphereTypes atmosphereType ):
        atmosphereType_( atmosphereType ){ }

    //! Destructor
    virtual ~AtmosphereSettings( ){ }

    //! Function to return type of atmosphere model that is to be created.
    /*!
     *  Function to return type of atmosphere model that is to be created.
     *  \return Type of atmosphere model that is to be created.
     */
    AtmosphereTypes getAtmosphereType( ){ return atmosphereType_; }

    //! Function to return settings for the atmosphere's wind model.
    /*!
     *  Function to return settings for the atmosphere's wind model.
     *  \return Settings for the atmosphere's wind model.
     */
    boost::shared_ptr< WindModelSettings > getWindSettings( )
    {
        return windSettings_;
    }

    //! Function to (re)set settings for the atmosphere's wind model.
    /*!
     *  Function to (re)set settings for the atmosphere's wind model.
     *  \param windSettings Settings for the atmosphere's wind model.
     */
    void setWindSettings( const boost::shared_ptr< WindModelSettings > windSettings )
    {
        windSettings_ = windSettings;
    }

private:

    //!  Type of atmosphere model that is to be created.
    AtmosphereTypes atmosphereType_;

    //! Settings for the atmosphere's wind model.
    boost::shared_ptr< WindModelSettings > windSettings_;
};

//! AtmosphereSettings for defining an exponential atmosphere.
class ExponentialAtmosphereSettings: public AtmosphereSettings
{
public:
    //! Constructor.
    /*!
     *  Constructor.
     *  \param densityScaleHeight Scale height for density profile of atmosphere.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param densityAtZeroAltitude Atmospheric density at ground level.
     *  \param specificGasConstant Specific gas constant for (constant) atmospheric chemical
     *  composition.
     *  \param ratioOfSpecificHeats Ratio of specific heats for (constant) atmospheric chemical
     *  composition.
     */

    ExponentialAtmosphereSettings(
            const double densityScaleHeight, const double constantTemperature,
            const double densityAtZeroAltitude,
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( exponential_atmosphere ),
        densityScaleHeight_( densityScaleHeight ), constantTemperature_( constantTemperature ),
        densityAtZeroAltitude_( densityAtZeroAltitude ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats )
    { }

    //! Function to return scale heigh for density profile of atmosphere.
    /*!
     *  Function to return scale heigh for density profile of atmosphere.
     *  \return Scale heigh for density profile of atmosphere.
     */
    double getDensityScaleHeight( ){ return densityScaleHeight_; }

    //! Function to return constant atmospheric temperature.
    /*!
     *  Function to return constant atmospheric temperature.
     *  \return Constant atmospheric temperature.
     */
    double getConstantTemperature( ){ return constantTemperature_; }

    //! Function to return atmospheric density at ground level.
    /*!
     *  Function to return atmospheric density at ground level.
     *  \return Atmospheric density at ground level.
     */
    double getDensityAtZeroAltitude( ){ return densityAtZeroAltitude_; }

    //! Function to return specific gas constant for (constant) atmospheric chemical
    /*!
     *  Function to return specific gas constant for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getSpecificGasConstant( ){ return specificGasConstant_; }

    //! Function to return ratio of specific heats for (constant) atmospheric chemical
    /*!
     *  Function to return ratio of specific heats for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getRatioOfSpecificHeats( ){ return ratioOfSpecificHeats_; }

private:

    //! Scale heigh for density profile of atmosphere.
    double densityScaleHeight_;

    //! Constant atmospheric temperature.
    double constantTemperature_;

    //! Atmospheric density at ground level.
    double densityAtZeroAltitude_;

    //! Specific gas constant for (constant) atmospheric chemical
    double specificGasConstant_;

    //! Ratio of specific heats for (constant) atmospheric chemical
    double ratioOfSpecificHeats_;
};


//! AtmosphereSettings for defining an NRLMSISE00 atmosphere reading space weather data from a text file.
class NRLMSISE00AtmosphereSettings: public AtmosphereSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param spaceWeatherFile File containing space weather data, as in
     *  https://celestrak.com/SpaceData/sw19571001.txt
     */
    NRLMSISE00AtmosphereSettings( const std::string& spaceWeatherFile ):
        AtmosphereSettings( nrlmsise00 ), spaceWeatherFile_( spaceWeatherFile ){ }

    //! Function to return file containing space weather data.
    /*!
     *  Function to return file containing space weather data.
     *  \return Filename containing space weather data.
     */
    std::string getSpaceWeatherFile( ){ return spaceWeatherFile_; }

private:

    //! File containing space weather data.
    /*!
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;
};


//! AtmosphereSettings for defining an atmosphere with tabulated data from file.
class TabulatedAtmosphereSettings: public AtmosphereSettings
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *          independent and dependent parameters needs to be specified in the independentVariablesNames and
     *          dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *          will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *          use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames =
    { altitude_dependent_atmosphere },
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames =
    { density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere },
                                 const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                 const double ratioOfSpecificHeats = 1.4,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling = { },
                                 const std::vector< double >& defaultExtrapolationValue = { } ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ),
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    { }

    //! Constructor with single boundary handling parameters.
    /*!
     *  Constructor with single boundary handling parameters. The specifier is assumed to be the same for
     *  each (in)dependent variable.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *          independent and dependent parameters needs to be specified in the independentVariablesNames and
     *          dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *          will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *          use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const double specificGasConstant,
                                 const double ratioOfSpecificHeats,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement< double >::getAdditionIdentity( ) ):
        TabulatedAtmosphereSettings( atmosphereTableFile, independentVariablesNames, dependentVariablesNames,
                                     specificGasConstant, ratioOfSpecificHeats,
                                     std::vector< interpolators::BoundaryInterpolationType >(
                                         independentVariablesNames.size( ), boundaryHandling ),
                                     std::vector< double >( dependentVariablesNames.size( ),
                                                            defaultExtrapolationValue ) ){ }

    //! Constructor compatible with old version.
    /*!
     *  Constructor compatible with old version.
     *  \param atmosphereTableFile File containing atmospheric properties.
     *          The file name of the atmosphere table. The file should contain four columns of data,
     *          containing altitude (first column), and the associated density, pressure and density values
     *          in the second, third and fourth columns.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     */
    TabulatedAtmosphereSettings( const std::string& atmosphereTableFile,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = {
                    density_dependent_atmosphere, pressure_dependent_atmosphere, temperature_dependent_atmosphere },
                                 const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                 const double ratioOfSpecificHeats = 1.4,
                                 const interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_boundary_value,
                                 const double defaultExtrapolationValue = IdentityElement< double >::getAdditionIdentity( ) ) :
        TabulatedAtmosphereSettings( { { 0, atmosphereTableFile } }, { altitude_dependent_atmosphere },
                                     dependentVariablesNames, specificGasConstant,
                                     ratioOfSpecificHeats, { boundaryHandling },
                                     std::vector< double >( dependentVariablesNames.size( ),
                                                            defaultExtrapolationValue ) ){ }

    //! Constructor with no specific gas constant nor ratio of specific heats.
    /*!
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector).
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *          independent and dependent parameters needs to be specified in the independentVariablesNames and
     *          dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *          will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *          use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                                 const std::vector< double >& defaultExtrapolationValue ) :
        TabulatedAtmosphereSettings( atmosphereTableFile, independentVariablesNames, dependentVariablesNames,
                                     physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4, boundaryHandling,
                                     defaultExtrapolationValue ){ }

    //! Constructor with no specific gas constant nor ratio of specific heats, and with
    //! single boundary handling parameters.
    /*!
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector). Only one boundary handling parameter is specified, which is then repeated for
     *  dimension.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *          independent and dependent parameters needs to be specified in the independentVariablesNames and
     *          dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *          will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *          use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement< double >::getAdditionIdentity( ) ) :
        TabulatedAtmosphereSettings( atmosphereTableFile, independentVariablesNames, dependentVariablesNames,
                                     physical_constants::SPECIFIC_GAS_CONSTANT_AIR, 1.4,
                                     std::vector< interpolators::BoundaryInterpolationType >(
                                         independentVariablesNames.size( ), boundaryHandling ),
                                     std::vector< double >( dependentVariablesNames.size( ),
                                                            defaultExtrapolationValue ) ){ }

    //! Function to return file containing atmospheric properties.
    /*!
     *  Function to return file containing atmospheric properties.
     *  \return Map of filenames containing atmospheric properties.
     */
    std::map< int, std::string > getAtmosphereFile( ){ return atmosphereFile_; }

    //! Function to return file containing atmospheric properties.
    /*!
     *  Function to return file containing atmospheric properties.
     *  \return Filename containing atmospheric properties.
     */
    std::string getAtmosphereFile( const unsigned int fileIndex ){ return atmosphereFile_.at( fileIndex ); }

    //! Function to return independent variables names.
    /*!
     *  Function to return independent variables names.
     *  \return Independent variables.
     */
    std::vector< AtmosphereIndependentVariables > getIndependentVariables( ){ return independentVariables_; }

    //! Function to return dependent variables names.
    /*!
     *  Function to return dependent variables names.
     *  \return Dependent variables.
     */
    std::vector< AtmosphereDependentVariables > getDependentVariables( ){ return dependentVariables_; }

    //! Function to return specific gas constant of the atmosphere.
    /*!
     *  Function to return specific gas constant of the atmosphere.
     *  \return Specific gas constant of the atmosphere.
     */
    double getSpecificGasConstant( ){ return specificGasConstant_; }

    //! Function to return ratio of specific heats of the atmosphere.
    /*!
     *  Function to return ratio of specific heats of the atmosphere.
     *  \return Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double getRatioOfSpecificHeats( ){ return ratioOfSpecificHeats_; }

    //! Function to return boundary handling method.
    /*!
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
     std::vector< interpolators::BoundaryInterpolationType > getBoundaryHandling( ){ return boundaryHandling_; }

    //! Function to return default extrapolation value.
    /*!
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< double > getDefaultExtrapolationValue( ){ return defaultExtrapolationValue_; }

private:

    //! File containing atmospheric properties.
    /*!
     *  File containing atmospheric properties, file should contain
     *  columns of atmospheric data with at least density, pressure and temperature,
     *  (whose order is specified in dependentVariables), and with at least one
     *  indendent variables.
     */
    std::map< int, std::string > atmosphereFile_;

    //! A vector of strings containing the names of the independent variables contained in the atmosphere file
    /*!
     * A vector of strings containing the names of the independent variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereIndependentVariables > independentVariables_;

    //! A vector of strings containing the names of the variables contained in the atmosphere file
    /*!
     * A vector of strings containing the names of the variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereDependentVariables > dependentVariables_;

    //! Specific gas constant of the atmosphere.
    /*!
     * Specific gas constant of the atmosphere.
     */
    double specificGasConstant_;

    //! Ratio of specific heats of the atmosphere at constant pressure and constant volume.
    /*!
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double ratioOfSpecificHeats_;

    //! Behavior of interpolator when independent variable is outside range.
    /*!
     *  Behavior of interpolator when independent variable is outside range.
     */
    std::vector< interpolators::BoundaryInterpolationType > boundaryHandling_;

    //! Default value to be used for extrapolation.
    /*!
     *  Default value to be used for extrapolation.
     */
    std::vector< double > defaultExtrapolationValue_;
};

//! Function to create a wind model.
/*!
 *  Function to create a wind model based on model-specific settings for the wind model.
 *  \param windSettings Settings for the wind model that is to be created, defined
 *  a pointer to an object of class (derived from) WindModelSettings.
 *  \param body Name of the body for which the wind model is to be created.
 *  \return Wind model created according to settings in windSettings.
 */
boost::shared_ptr< WindModel > createWindModel(
        const boost::shared_ptr< WindModelSettings > windSettings,
        const std::string& body );

//! Function to create an atmosphere model.
/*!
 *  Function to create an atmosphere model based on model-specific settings for the atmosphere.
 *  \param atmosphereSettings Settings for the atmosphere model that is to be created, defined
 *  a pointer to an object of class (derived from) AtmosphereSettings.
 *  \param body Name of the body for which the atmosphere model is to be created.
 *  \return Atmosphere model created according to settings in atmosphereSettings.
 */
boost::shared_ptr< AtmosphereModel > createAtmosphereModel(
        const boost::shared_ptr< AtmosphereSettings > atmosphereSettings,
        const std::string& body );

} // namespace simulation_setup

} // namespace tudat


#endif // TUDAT_CREATEATMOSPHEREMODEL_H
