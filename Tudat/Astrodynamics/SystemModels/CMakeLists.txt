 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the header files.
set(SYSTEMMODELS_HEADERS 
  "${SRCROOT}${SYSTEMMODELSDIR}/controlSystem.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/engineModel.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/extraFunctions.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/guidanceSystem.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/instrumentsModel.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/navigationSystem.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/onboardComputerModel.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/vehicleSystems.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/UnitTests/generators.h"
)

set(SYSTEMMODELS_SOURCES
  "${SRCROOT}${SYSTEMMODELSDIR}/controlSystem.cpp"
  "${SRCROOT}${SYSTEMMODELSDIR}/extraFunctions.cpp"
  "${SRCROOT}${SYSTEMMODELSDIR}/guidanceSystem.cpp"
  "${SRCROOT}${SYSTEMMODELSDIR}/instrumentsModel.cpp"
  "${SRCROOT}${SYSTEMMODELSDIR}/navigationSystem.cpp"
)

# Add static libraries.
add_library(tudat_system_models STATIC ${SYSTEMMODELS_SOURCES} ${SYSTEMMODELS_HEADERS})
setup_tudat_library_target(tudat_system_models "${SRCROOT}${SYSTEMMODELSDIR}")
set_target_properties(tudat_system_models PROPERTIES LINKER_LANGUAGE CXX)

# Add unit tests
add_executable(test_InstrumentsModel "${SRCROOT}${SYSTEMMODELSDIR}/UnitTests/unitTestInstrumentsModel.cpp")
setup_custom_test_program(test_InstrumentsModel "${SRCROOT}${SYSTEMMODELSDIR}")
target_link_libraries(test_InstrumentsModel ${TUDAT_PROPAGATION_LIBRARIES} ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_NavigationSystem "${SRCROOT}${SYSTEMMODELSDIR}/UnitTests/unitTestNavigationSystem.cpp")
setup_custom_test_program(test_NavigationSystem "${SRCROOT}${SYSTEMMODELSDIR}")
target_link_libraries(test_NavigationSystem ${TUDAT_PROPAGATION_LIBRARIES} ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_GuidanceSystem "${SRCROOT}${SYSTEMMODELSDIR}/UnitTests/unitTestGuidanceSystem.cpp")
setup_custom_test_program(test_GuidanceSystem "${SRCROOT}${SYSTEMMODELSDIR}")
target_link_libraries(test_GuidanceSystem ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_ControlSystem "${SRCROOT}${SYSTEMMODELSDIR}/UnitTests/unitTestControlSystem.cpp")
setup_custom_test_program(test_ControlSystem "${SRCROOT}${SYSTEMMODELSDIR}")
target_link_libraries(test_ControlSystem ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
