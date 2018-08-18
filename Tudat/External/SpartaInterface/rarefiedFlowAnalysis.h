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
 *      Klothakis, A. and Nikolos, I., “Modeling of Rarefied Hypersonic Flows Using the Massively
 *        Parallel DSMC Kernel “SPARTA”,” in 8th GRACM International Congress on Computational Mechanics,
 *        Volos, Greece, July 2015.
 *      Plimpton, S. and Gallis, M., SPARTA Users Manual, Sandia National Laboratories, United States
 *        Department of Energy, July 2017.
 *      Liechty, D., “Aeroheating Analysis for the Mars Reconnaissance Orbiter with Comparison to Flight Data,”
 *        Journal of Spacecraft and Rockets, vol. 44, no. 6, pp. 1226–1231, 2007.
 */

#ifndef TUDAT_RAREFIED_FLOW_ANALYSIS_H
#define TUDAT_RAREFIED_FLOW_ANALYSIS_H

#include <map>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/External/SpartaInterface/spartaInterface.h"

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientGenerator.h"
#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace aerodynamics
{

//! Returns default values of altitude for use in RarefiedFlowAnalysis.
/*!
 *  Returns default values of altitude for use in RarefiedFlowAnalysis.
 */
std::vector< double > getDefaultRarefiedFlowAltitudePoints(
        const std::string& targetPlanet = "Earth" );

//! Returns default values of Mach number for use in RarefiedFlowAnalysis.
/*!
 *  Returns default values of Mach number for use in RarefiedFlowAnalysis.
 */
std::vector< double > getDefaultRarefiedFlowMachPoints(
        const std::string& machRegime = "Full" );

//! Returns default values of angle of attack for use in RarefiedFlowAnalysis.
/*!
 *  Returns default values of angle of attack for use in RarefiedFlowAnalysis.
 */
std::vector< double > getDefaultRarefiedFlowAngleOfAttackPoints(
        const std::string& angleOfAttackRegime = "Reduced" );

//! Class for aerodynamic analysis in rarefied flow using the SPARTA DSMC method.
/*!
 *  Class for aerodynamic analysis in rarefied flow using the SPARTA DSMC (Direct Simulation
 *  Monte Carlo) method. This method uses a Monte Carlo simulation, to determine the pressure and
 *  shear forces acting on each element of the vehicle. These values are output by default every 200
 *  time steps, and are used to compute the average pressure and friction coefficients on each surface
 *  element, which are then translated to aerodynamic coefficients for the whole surface. One can
 *  find a description of the SPARTA software in references [Klothakis and Nikolos, 2015, Plimpton and
 *  Gallis, 2017], where the second reference is the official user manual. The user should also pay careful
 *  attention to the requirements for the geometry of the vehicle to be analyzed.
 */
class RarefiedFlowAnalysis: public AerodynamicCoefficientGenerator< 3, 6 >, public sparta_interface::SpartaInterface
{
public:

    // Using statements
    using AerodynamicCoefficientGenerator< 3, 6 >::dataPointsOfIndependentVariables_;

    //! Default constructor.
    /*!
     *  Default constructor for SPARTA rarefied flow simulation.
     *  \param SpartaExecutable Path to executable for SPARTA simulation.
     *  \param dataPointsOfIndependentVariables Vector of vectors, with each subvector containing
     *          the data points of each of the independent variables for the coefficient generation.
     *          The physical meaning of each of the three independent variables is: 0 = altitude,
     *          1 = Mach number, 1 = angle of attack. Each of the subvectors must be sorted in ascending order.
     *  \param simulationGases String of gases making up the atmosphere of the planet, to be used for
     *          the simulation.
     *  \param atmosphereModel Pointer to the atmosphere model of the planet. Should provide information on
     *          density, pressure, temperature, gas constant and specific heat ratio.
     *  \param geometryFileUser Path to the file describing the geometry of the vehicle, where
     *          the surface elements are discretized as triangles.
     *  \param referenceArea Reference area used to non-dimensionalize aerodynamic forces
     *          and moments.
     *  \param referenceLength Reference length used to non-dimensionalize aerodynamic moments.
     *  \param referenceAxis Index of main axis of the vehicle (i.e., axis opposite in direction to
     *          incoming flow, when angles of attack and sideslip are zero).
     *  \param momentReferencePoint Reference point wrt which aerodynamic moments are calculated.
     *  \param gridSpacing Grid size for simulation environment, used to define the size of each cell, and
     *          the number of cells in the environment.
     *  \param simulatedParticlesPerCell Number of simulated particles per cell.
     *  \param wallTemperature Temperature of the surface of the vehicle (default value is 300 K [Liechty, 2007]).
     *  \param accommodationCoefficient Accommodation coefficient of the surface of the vehicle. This
     *          value indicates the degree of diffusivity during molecular-surface collisions (default value
     *          is 1.0, i.e., diffuse reflection).
     *  \param printProgressInCommandWindow Boolean to toggle appearance of SPARTA output in command
     *          window (default is false).
     *  \param MpiExecutable Path to executable for multi-processor computing software (default is none).
     *  \param numberOfCores Number of cores to be allocated for multi-processor computing of
     *          SPARTA simulation (default is 0). Note that this value can only be used if the path to MPI
     *          is set. Furthermore, the number of cores has to be an integer larger or equal to one.
     */
    RarefiedFlowAnalysis(
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            boost::shared_ptr< TabulatedAtmosphere > atmosphereModel,
            const std::string& simulationGases,
            const std::string& geometryFileUser,
            const double referenceArea,
            const double referenceLength,
            const int referenceAxis,
            const Eigen::Vector3d& momentReferencePoint,
            const double gridSpacing,
            const double simulatedParticlesPerCell,
            const double wallTemperature = 300.0,
            const double accommodationCoefficient = 1.0,
            const bool printProgressInCommandWindow = false,
            const std::string& SpartaExecutable = "~/sparta/src/spa_mpi",
            const std::string& MpiExecutable = "mpirun",
            const unsigned int numberOfCores = 14 );

    //! Default destructor.
    /*!
     *  Default destructor.
     */
    ~RarefiedFlowAnalysis( ) { }

    //! Get aerodynamic coefficients.
    /*!
     *  Returns aerodynamic coefficients.
     *  The physical meaning of each of the three independent variables is: 0 = altitude,
     *  1 = Mach number, 1 = angle of attack.
     *  \param independentVariables Array of values of independent variable
     *          indices in dataPointsOfIndependentVariables_.
     *  \return vector of coefficients at specified independent variable indices.
     */
    Eigen::Vector6d getAerodynamicCoefficientsDataPoint(
            const boost::array< int, 3 > independentVariables );

private:

    //! Generate aerodynamic database.
    /*!
     *  Generates aerodynamic database, by running the SPARTA simulation via command line (standard
     *  library function std::system), and reads output of simulation to extract values of aerodynamic
     *  coefficients, as a function of altitude, Mach number and angle of attack.
     */
    void generateCoefficients( );

};

//! Typedef for shared-pointer to RarefiedFlowAnalysis object.
typedef boost::shared_ptr< RarefiedFlowAnalysis > RarefiedFlowAnalysisPointer;

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_RAREFIED_FLOW_ANALYSIS_H
