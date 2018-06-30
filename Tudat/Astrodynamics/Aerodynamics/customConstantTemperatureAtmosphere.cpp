#include "Tudat/Astrodynamics/Aerodynamics/customConstantTemperatureAtmosphere.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

namespace tudat
{

namespace aerodynamics
{

//! First atmosphere model, based on exponential atmosphere.
double exponentialAtmosphereModel( const double altitude, const double longitude, const double latitude, const double time,
                                   const double densityAtReferenceAltitude, const double referenceAltitude, const double scaleHeight )
{
    // Compute density
    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );
    return densityAtReferenceAltitude * std::exp( ( referenceAltitude - altitude ) / scaleHeight );
}

//! Second atmosphere model, based on a three longitudinal waves model.
double threeWaveAtmosphereModel( const double altitude, const double longitude, const double latitude, const double time,
                                 const double densityAtReferenceAltitude, const double referenceAltitude, const double scaleHeight,
                                 const double uncertaintyFactor, const double dustStormFactor )
{
    // Compute density
    double waveModelVariation = uncertaintyFactor + dustStormFactor +
            0.1 * std::sin( longitude ) + // first longitudinal wave
            0.2 * std::sin( longitude - unit_conversions::convertDegreesToRadians( 50.0 ) ) + // second longitudinal wave
            0.1 * std::sin( longitude - unit_conversions::convertDegreesToRadians( 55.0 ) ); // third longitudinal wave
    return exponentialAtmosphereModel( altitude, longitude, latitude, time, densityAtReferenceAltitude,
                                       referenceAltitude, scaleHeight ) * waveModelVariation;
}

//! Third atmosphere model, based on three constant scale height atmospheres.
double threeTermAtmosphereModel( const double altitude, const double longitude, const double latitude, const double time,
                                 const double densityAtReferenceAltitude, const double referenceAltitude, const double scaleHeight,
                                 const std::vector< double >& modelWeights )
{
    using namespace mathematical_constants;

    // Compute density
    TUDAT_UNUSED_PARAMETER( longitude );
    TUDAT_UNUSED_PARAMETER( latitude );
    TUDAT_UNUSED_PARAMETER( time );
    return densityAtReferenceAltitude *
            std::exp( modelWeights.at( 0 ) * ( altitude - referenceAltitude ) / scaleHeight + // first CSH model
                      modelWeights.at( 1 ) * std::cos( 2.0 * PI * ( altitude - referenceAltitude ) / scaleHeight ) + // second CSH model
                      modelWeights.at( 2 ) * std::sin( 2.0 * PI * ( altitude - referenceAltitude ) / scaleHeight ) ); // third CSH model
}

//! Constructor which uses one of the built-in density functions as input.
CustomConstantTemperatureAtmosphere::CustomConstantTemperatureAtmosphere(
        const AvailableConstantTemperatureAtmosphereModels densityFunctionIdentifier,
        const double constantTemperature,
        const double specificGasConstant,
        const double ratioOfSpecificHeats,
        const std::vector< double >& modelSpecificParameters ) :
    constantTemperature_( constantTemperature ),
    specificGasConstant_( specificGasConstant ),
    ratioOfSpecificHeats_( ratioOfSpecificHeats )
{
    // Set density function based on user-provided data
    switch ( densityFunctionIdentifier )
    {
    case exponential_atmosphere_model:
    {
        // Check that the number of input parameters is correct
        if ( modelSpecificParameters.size( ) != 3 )
        {
            throw std::runtime_error( "Error while creating custom constant temperature atmosphere model. The "
                                      "number of input model-dependent parameters is incorrect. Number of "
                                      "parameters given " + modelSpecificParameters.size( ) + ". Number required: 3." );
        }

        // Set density function
        densityFunction_ = boost::bind( &exponentialAtmosphereModel, _1, _2, _3, _4,
                                        modelSpecificParameters.at( 0 ), modelSpecificParameters.at( 1 ),
                                        modelSpecificParameters.at( 2 ) );
        break;
    }
    case three_wave_atmosphere_model:
    {
        // Check that the number of input parameters is correct
        if ( modelSpecificParameters.size( ) != 5 )
        {
            throw std::runtime_error( "Error while creating custom constant temperature atmosphere model. The "
                                      "number of input model-dependent parameters is incorrect. Number of "
                                      "parameters given " + modelSpecificParameters.size( ) + ". Number required: 5." );
        }

        // Set density function
        densityFunction_ = boost::bind( &threeWaveAtmosphereModel, _1, _2, _3, _4,
                                        modelSpecificParameters.at( 0 ), modelSpecificParameters.at( 1 ),
                                        modelSpecificParameters.at( 2 ), modelSpecificParameters.at( 3 ),
                                        modelSpecificParameters.at( 4 ) );
        break;
    }
    case three_term_atmosphere_model:
    {
        // Check that the number of input parameters is correct
        if ( modelSpecificParameters.size( ) != 4 )
        {
            throw std::runtime_error( "Error while creating custom constant temperature atmosphere model. The "
                                      "number of input model-dependent parameters is incorrect. Number of "
                                      "parameters given " + modelSpecificParameters.size( ) + ". Number required: 4." );
        }
        if ( modelSpecificParameters.at( 3 ).size( ) != 3 )
        {
            throw std::runtime_error( "Error while creating custom constant temperature atmosphere model. The "
                                      "number of input weights for the three-term model is incorrect. Number of "
                                      "weights given " + modelSpecificParameters.at( 3 ).size( ) + ". Number required: 3." );
        }

        // Set density function
        densityFunction_ = boost::bind( &threeTermAtmosphereModel, _1, _2, _3, _4,
                                        modelSpecificParameters.at( 0 ), modelSpecificParameters.at( 1 ),
                                        modelSpecificParameters.at( 2 ), modelSpecificParameters.at( 3 ) );
        break;
    }
    default:
        throw std::runtime_error( "Error while creating custom constant temperature atmosphere model. The "
                                  "model selected is not recognized." );
    }
}

} // namespace aerodynamics

} // namespace tudat
