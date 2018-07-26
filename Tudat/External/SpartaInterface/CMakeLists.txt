 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(SPARTAINTERFACE_SOURCES
  "${SRCROOT}${EXTERNALDIR}/SpartaInterface/spartaInterface.cpp"
  "${SRCROOT}${EXTERNALDIR}/SpartaInterface/rarefiedFlowAnalysis.cpp"
)

# Set the header files.
set(SPARTAINTERFACE_HEADERS
  "${SRCROOT}${EXTERNALDIR}/SpartaInterface/spartaInterface.h"
  "${SRCROOT}${EXTERNALDIR}/SpartaInterface/rarefiedFlowAnalysis.h"
  "${SRCROOT}${EXTERNALDIR}/SpartaInterface/spartaInputOutput.h"
  "${SRCROOT}${EXTERNALDIR}/SpartaInterface/spartaDataReader.h"
)

# Add static libraries.
add_library(tudat_sparta_interface STATIC ${SPARTAINTERFACE_SOURCES} ${SPARTAINTERFACE_HEADERS})
setup_tudat_library_target(tudat_sparta_interface "${SRCROOT}${SPARTAINTERFACEDIR}")

# Add unit tests.
