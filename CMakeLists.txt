##
#  CMake script for the waveguide-with-refractor project
#  (currently: mesh test only)
##

SET(TARGET "main")                 # solver target, not needed yet
#SET(TARGET "test_mesh")

SET(TARGET_SRC
  #tests/mesh_test/mesh_test.cc       # adjust if your test file is elsewhere
  src/mesh_builder.cc
  src/main.cc                          
  src/parameters.cc  
  src/post_processor.cc   
  src/wave_guide_problem.cc                  
  )

CMAKE_MINIMUM_REQUIRED(VERSION 3.13.4)

FIND_PACKAGE(deal.II 9.3.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})

# Gmsh (C++ API)
FIND_LIBRARY(GMSH_LIB gmsh REQUIRED)
FIND_PATH(GMSH_INC gmsh.h REQUIRED)

# header folder and Gmsh headers
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/includes ${GMSH_INC})
LINK_LIBRARIES(${GMSH_LIB})

DEAL_II_INVOKE_AUTOPILOT()