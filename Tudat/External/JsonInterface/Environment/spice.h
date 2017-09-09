/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_SPICE_H
#define TUDAT_JSONINTERFACE_SPICE_H

#include <Tudat/External/SpiceInterface/spiceInterface.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace json_interface
{

/*
// SimulationType

//! Frequently-used simulations.
enum SimulationType
{
    customSimulation,
    singlePerturbedBody
};

//! Map of `SimulationType` string representations.
static std::map< SimulationType, std::string > simulationTypes =
{
    { customSimulation, "custom" },
    { singlePerturbedBody, "singlePerturbedBody" }
};

//! `SimulationType` not supported by `json_interface`.
static std::vector< SimulationType > unsupportedSimulationTypes = {  };

//! Convert `SimulationType` to `json`.
inline void to_json( json& jsonObject, const SimulationType& simulationType )
{
    jsonObject = json_interface::stringFromEnum( simulationType, simulationTypes );
}

//! Convert `json` to `SimulationType`.
inline void from_json( const json& jsonObject, SimulationType& simulationType )
{
    simulationType = json_interface::enumFromString( jsonObject, simulationTypes );
}
*/


// SpiceSettings

/*
//! Get the set of spice kernels to be used for a SimulationType.
std::vector< boost::filesystem::path > getSpiceKernels( const SimulationType simulationType );
*/

//! Class containing the settings for Spice used in a simulation.
/*!
 * Class containing the settings for Spice used in a simulation.
 */
class SpiceSettings
{
public:
    /*
    //! Constructor.
    SpiceSettings( const SimulationType simulationType ) :
        kernels_( getSpiceKernels( simulationType ) ) { }
    */

    //! Constructor with a vector of Spice kernels to be used.
    /*!
     * @copybrief SpiceSettings
     * \param kernels
     */
    SpiceSettings( const std::vector< boost::filesystem::path >& kernels ) : kernels_( kernels ) { }

    //! Destructor.
    virtual ~SpiceSettings( ) { }


    //! Vector containing the paths to the spice kernel files to be loaded.
    std::vector< boost::filesystem::path > kernels_;

    //! Whether all the data from the Spice kernels should be preloaded before the simulation for the interval start
    //! epoch to end epoch (true), or whether the data from Spice should be accessed on request at every step (false).
    /*!
     * Whether all the data from the Spice kernels should be preloaded before the simulation for the interval start
     * epoch to end epoch (true), or whether the data from Spice should be accessed on request at every step (false).
     * <br/>
     * Preloading Spice data generally results in faster propagations, unless:
     * <br/>
     * <ul>
     *  <li>The simulation ends much earlier than the specified maximum simulation end epoch.</li>
     *  <li>The integrator step-size is very large (in the order of several hours or days).</li>
     * </ul>
     */
    bool preloadKernels_ = true;

    //! Offsets for the interval for which the spice kernels are to be preloaded.
    /*!
     * Offsets for the interval for which the spice kernels are to be preloaded.
     * <br/>
     * The kernels will be preloaded for the interval:
     * `[ initialEpoch + preloadOffsets_.first, finalEpoch + preloadOffsets_.second ]`
     * \remark Ignored if SpiceSettings::preloadKernels_ is set to `false`.
     * \remark If not specified, the used values are 10 * interpolationStep_.
     */
    std::pair< double, double > preloadOffsets_ = { TUDAT_NAN, TUDAT_NAN };

    //! Step-size for the interpolated Spice ephemeris.
    /*!
     * Step-size for the interpolated Spice ephemeris. Ignored if preloadKernels_ set to false.
     */
    double interpolationStep_ = 300.0;

    //! Get initial offset for the interpolated Spice ephemeris.
    /*!
     * @copybrief getInitialOffset
     * \remark If not defined by the user (i.e. is NaN), returns `-10 * interpolationStep_`.
     * \return Initial offset for the interpolated Spice ephemeris.
     */
    double getInitialOffset( )
    {
        if ( isNaN( preloadOffsets_.first ) )
        {
            return -10.0 * interpolationStep_;
        }
        else
        {
            return preloadOffsets_.first;
        }
    }

    //! Get final offset for the interpolated Spice ephemeris.
    /*!
     * @copybrief getFinalOffset
     * \remark If not defined by the user (i.e. is NaN), returns `10 * interpolationStep_`.
     * \return Final offset for the interpolated Spice ephemeris.
     */
    double getFinalOffset( )
    {
        if ( isNaN( preloadOffsets_.second ) )
        {
            return 10.0 * interpolationStep_;
        }
        else
        {
            return preloadOffsets_.second;
        }
    }
};

//! Create a `json` object from a shared pointer to a `SpiceSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< SpiceSettings >& spiceSettings );

//! Create a shared pointer to a `SpiceSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< SpiceSettings >& spiceSettings );


//! Load in Tudat the Spice kernels specified in \p spiceSettings.
/*!
 * @copybrief loadSpiceKernels
 * \remark Clears any Spice kernel loaded previously.
 * \remark If \p spiceSettings is `NULL`, no kernels are loaded.
 * \param spiceSettings The Spice settings containing the paths to the kernels to be loaded.
 */
void loadSpiceKernels( const boost::shared_ptr< SpiceSettings >& spiceSettings );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SPICE_H
