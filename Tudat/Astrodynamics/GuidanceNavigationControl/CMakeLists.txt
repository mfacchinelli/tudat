 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(GUIDANCENAVIGATIONCONTROL_SOURCES
  "${SRCROOT}${GUIDANCENAVIGATIONCONTROLDIR}/controlSystem.cpp"
  "${SRCROOT}${GUIDANCENAVIGATIONCONTROLDIR}/guidanceSystem.cpp"
  "${SRCROOT}${GUIDANCENAVIGATIONCONTROLDIR}/navigationSystem.cpp"
)

# Set the header files.
set(GUIDANCENAVIGATIONCONTROL_HEADERS
  "${SRCROOT}${GUIDANCENAVIGATIONCONTROLDIR}/controlSystem.h"
  "${SRCROOT}${GUIDANCENAVIGATIONCONTROLDIR}/guidanceSystem.h"
  "${SRCROOT}${GUIDANCENAVIGATIONCONTROLDIR}/navigationSystem.h"
)

# Add static libraries.
add_library(tudat_michele_gnc STATIC ${GUIDANCENAVIGATIONCONTROL_SOURCES} ${GUIDANCENAVIGATIONCONTROL_HEADERS})
setup_tudat_library_target(tudat_michele_gnc "${SRCROOT}{GUIDANCENAVIGATIONCONTROLDIR}")
