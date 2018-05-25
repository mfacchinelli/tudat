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

#ifndef TUDAT_SPARTA_INTERFACE_H
#define TUDAT_SPARTA_INTERFACE_H

#include <map>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/External/SpartaInterface/rarefiedFlowAnalysis.h"
#include "Tudat/External/SpartaInterface/spartaDataReader.h"
#include "Tudat/External/SpartaInterface/spartaInputOutput.h"

namespace tudat
{

namespace sparta_interface
{

//! Enumaration of key values for map of atmospheric conditions.
enum AtmosphericConditionVariables
{
    density_index = 0,
    pressure_index = 1,
    temperature_index = 2,
    speed_of_sound_index = 3,
    number_density_index = 4
};

//! Function to sort the rows of a matrix, based on the specified column and specified order.
/*!
 *  Function to sort the rows of a matrix, based on the specified column and specified order. Note that
 *  this function only works if the column to be sorted is made up of consecutive integers (in a scrambled
 *  order, of course).
 *  \param matrixToBeSorted Matrix that has to be sorted.
 *  \param referenceColumn Column to be used as reference for the sorting.
 *  \param descendingOrder Boolean to toggle sorting in descending order.
 *  \return Matrix where the rows have been sorted based on the specified column number.
 */
Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > sortMatrixRows(
        const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >& matrixToBeSorted,
        const int referenceColumn, const bool descendingOrder = false );

//!
/*!
 *
 */
class SpartaInterface
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     */
    SpartaInterface( const std::string& SPARTAExecutable = "~/sparta/src/spa_mpi",
                     const std::string& MPIExecutable = "mpirun" ) :
        SPARTAExecutable_( SPARTAExecutable ), MPIExecutable_( MPIExecutable )
    {
        // Check executables
        checkExecutableValidity( );

        // Read SPARTA input template
        inputTemplate_ = input_output::readSpartaInputFileTemplate( );
    }

protected:

    //! Open and read geometry file for SPARTA simulation.
    /*!
     *  Open and read geometry file for SPARTA simulation, to extract information on the
     *  vehicle and surface elements dimensions and properties.
     *  \param geometryFileUser Path to the file describing the geometry of the vehicle, where
     *          the surface elements are discretized as triangles.
     *  \param momentReferencePoint Reference point wrt which aerodynamic moments are calculated.
     */
    void analyzeGeometryFile( const std::string& geometryFileUser, const double referenceArea,
                              const Eigen::Vector3d& momentReferencePoint );

    //! Retrieve simulation conditions based on input and geometry.
    /*!
     *  Retrieve simulation conditions based on input and geometry, including simulation environment boundaries,
     *  velocity of incoming flow, time step of simulation (set as 10 % of the time needed for a particle to
     *  traverse the simulation environment), and ratio of real-to-simulated particles.
     */
    void getSimulationConditions( );

    //! Run SPARTA simulation.
    /*!
     *  Run SPARTA simulation with conditions and input files specified by the user.
     */
    void runSpartaSimulation( const unsigned int h, const unsigned int m, const unsigned int a );

    //! Process SPARTA output.
    /*!
     *  Process SPARTA output by computing force and moment coefficients based on the pressure and shear force
     *  distribution over the vehicle.
     */
    Eigen::Vector6d processSpartaOutput( const double referenceArea, const double referenceLength,
                                         const unsigned int h, const unsigned int m );

    //! String of gases making up the atmosphere of the target planet.
    std::string simulationGases_;

    //! Reference axis for the aerodynamic analysis.
    /*!
     *  Reference axis for the aerodynamic analysis. This axis is used to set very important paramters for
     *  the simulation, e.g., the velocity direction, reference surface area, etc. It is thus fundamental
     *  to get this value right. The axis should be an integer (i.e., a signed integer), such that the flow of
     *  particles comes from the opposite direction of this axis.
     */
    int referenceAxis_;

    //! Grid size for simulation environment.
    /*!
     *  Grid size for simulation environment. Used to define the size of each cell, and the number of cells in
     *  the environment.
     */
    double gridSpacing_;

    //! Number of simulated particles per cell.
    double simulatedParticlesPerCell_;

    //! Temperature of surface of vehicle.
    double wallTemperature_;

    //! Accommodation coefficient of surface of vehicle.
    double accommodationCoefficient_;

    //! Boolean to toggle showing of progress in command window.
    /*!
     *  Boolean to toggle showing of progress in command window. SPARTA outputs information on the geometry and other environment
     *  details, and during the simulation it prints statistics on the progress. This value is set to false by default.
     */
    bool printProgressInCommandWindow_;

    //! Path to SPARTA executable.
    /*!
     *  Path to SPARTA executable. Note that SPARTA is an external software and needs to be compiled before
     *  it can be used in Tudat. See the instructions in the manual [2].
     */
    std::string SPARTAExecutable_;

    //! Path to open MPI executable.
    /*!
     *  Path to open MPI executable. Note that open MPI is an external software and needs to be compiled before
     *  it can be used in Tudat. See the instructions on the website https://www.open-mpi.org.
     */
    std::string MPIExecutable_;

    //! Number of cores to be used to run the simulation with open MPI.
    /*!
     *  Number of cores to be used to run the simulation with open MPI. Note that this can only be used if a path to the MPI
     *  executable has been set.
     */
    unsigned int numberOfCores_;

    //! Reference dimension for the aerodynamic analysis.
    /*!
     *  Reference dimension for the aerodynamic analysis. This value is based on the referenceAxis_ variable, but
     *  without the information on the sign, such that it can be used to access vector entries.
     */
    unsigned int referenceDimension_;

    //! List of points making up the vehicle geometry.
    /*!
     *  List of points making up the vehicle geometry. Column size is 3, since only triangular mesh are supported
     *  by SPARTA.
     */
    Eigen::Matrix< double, Eigen::Dynamic, 3 > shapePoints_;

    //! List of triangle vertices making up the vehicle geometry.
    /*!
     *  List of triangle vertices making up the vehicle geometry. Each element refers to a point from the shapePoints_
     *  list. Column size is 3, since only triangular mesh are supported by SPARTA.
     */
    Eigen::Matrix< int, Eigen::Dynamic, 3 > shapeTriangles_;

    //! Total number of points making up the vehicle geometry.
    int numberOfPoints_;

    //! Total number of triangles making up the vehicle geometry.
    int numberOfTriangles_;

    //! Maximum dimensions of the vehicle in x, y and z.
    Eigen::Vector3d maximumDimensions_;

    //! Minimum dimensions of the vehicle in x, y and z.
    Eigen::Vector3d minimumDimensions_;

    //! List of normal vector to each surface element.
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementSurfaceNormal_;

    //! List of surface area of each surface element.
    Eigen::Matrix< double, 1, Eigen::Dynamic > elementSurfaceArea_;

    //! List of moment arm for each surface element.
    /*!
     *  List of moment arm for each surface element. Note that the moment arm is defined as the distance
     *  of the centroid of each surface element, to a reference point, i.e., momentReferencePoint_.
     */
    Eigen::Matrix< double, 3, Eigen::Dynamic > elementMomentArm_;

    //! List of cross-sectional area of vehicle in x, y and z.
    Eigen::Vector3d shapeCrossSectionalArea_;

    //! Atmospheric conditions at altitudes specified by user in dataPointsOfIndependentVariables_.
    /*!
     *  Atmospheric conditions at altitudes specified by user in dataPointsOfIndependentVariables_. The key values
     *  are given by the enumeration AtmosphericConditionVariables, at the beginning of this file.
     */
    std::map< AtmosphericConditionVariables, std::vector< double > > atmosphericConditions_;

    //! List defining the simulation environment boundaries.
    /*!
     *  List defining the simulation environment boundaries. They are computed from the maximum and minimum
     *  dimensions of the vehicle, i.e., maximumDimensions_ and minimumDimensions_, by adding extra space in each
     *  direction (50 % of the vehicle size).
     */
    Eigen::Vector6d simulationBoundaries_;

    //! List of grid spacing for each dimension.
    /*!
     *  List of grid spacing for each dimension. They are determined by dividing the simulation dimensions (given
     *  by simulationBoundaries_) by the grid spacing, i.e., gridSpacing_.
     */
    Eigen::Vector3d simulationGrid_;

    //! List of freestream velocities for each altitude and Mach number.
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > freeStreamVelocities_;

    //! List of time steps for the simulation for each altitude and Mach number.
    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > simulationTimeStep_;

    //! List of ratios of real-to-simulated particles for each altitude.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > ratioOfRealToSimulatedParticles_;

    //! String containing the template of the input file for the SPARTA simulation.
    /*!
     *  String containing the template of the input file for the SPARTA simulation. This template contains
     *  format specifiers (such as %.3f) where a simulation variable is to be placed. This is done by filling
     *  the values input by the user and determined in previous parts of the code.
     */
    std::string inputTemplate_;

private:

    //! Function to check whether the input executables are indeed present in the system.
    /*!
     *  Function to check whether the input executables are indeed present in the system. If the SPARTA
     *  executable cannot be found, it will be created by downloading, installing and building the
     *  SPARTA source code. If the MPI executable cannot be found, the simulation will be run with only
     *  one core.
     */
    void checkExecutableValidity( )
    {
        // Check that SPARTA executable exists and create it otherwise
        if ( !boost::filesystem::exists( boost::filesystem::system_complete( SPARTAExecutable_ ) ) )
        {
            std::cerr << "SPARTA executable not found. "
                         "Cloning and building SPARTA from scratch." << std::endl;
            std::string cloneAndBuildSPARTACommandString =
                    "mkdir ~/sparta/; "
                    "cd ~/sparta/; git clone https://github.com/mfacchinelli/sparta.git; "
                    "cd ~/sparta/src/; make -j 14 mpi";
//            std::system( cloneAndBuildSPARTACommandString.c_str( ) );
        }

        // Check try running a dummy MPI example
        std::string testMPICommandString = "info " + MPIExecutable_;
        int systemStatus = std::system( testMPICommandString.c_str( ) );
        if ( systemStatus != 0 )
        {
            std::cerr <<  "Error in SPARTA interface. MPI executable not found. "
                          "Simulation will be run with one core only.";
            MPIExecutable_ = "";
        }

        // Create directory for output
        if ( !boost::filesystem::exists( input_output::getSpartaOutputPath( ) ) )
        {
            boost::filesystem::create_directories( input_output::getSpartaOutputPath( ) );
        }
    }

    //! Integer output by system command, specifying whether command was successfully executed.
    int systemStatus_;

    //! Velocity vector of air flow during simulation.
    Eigen::Vector3d velocityVector_;

    //! Path to files output by SPARTA simulation.
    std::string temporaryOutputFile_ = input_output::getSpartaOutputPath( ) + "/coeff";

    //! Extensions of files output by SPARTA simulation.
    std::vector< std::string > outputFileExtensions_ = { ".400", ".600", ".800", ".1000" };

    //! Matrix where output of SPARTA simulation is stored.
    Eigen::Matrix< double, Eigen::Dynamic, 7, Eigen::RowMajor > outputMatrix_;

    //! Matrix where mean pressure force vectors for each surface element are stored.
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanPressureValues_;

    //! Matrix where mean shear force vectors for each surface element are stored.
    Eigen::Matrix< double, 3, Eigen::Dynamic > meanShearValues_;

};

} // namespace sparta_interface

} // namespace tudat

#endif // TUDAT_SPARTA_INTERFACE_H
