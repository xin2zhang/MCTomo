cmake_minimum_required (VERSION 2.6)
project (examples)

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

# Select flags.
SET(CMAKE_CXX_FLAGS "-std=c++11")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -fno-finite-math-only")
SET(CMAKE_CXX_FLAGS_DEBUG  "-Wall -Wno-unused-local-typedefs -g")
if(APPLE)
    include_directories(AFTER "/opt/X11/include")
    link_directories( "/opt/X11/lib")
endif()

# Self-made includes
include_directories (${EXAMPLES_SOURCE_DIR}/../console)
include_directories (${EXAMPLES_SOURCE_DIR}/../ndgridmap)
include_directories (${EXAMPLES_SOURCE_DIR}/../fmm)
include_directories (${EXAMPLES_SOURCE_DIR}/../fm2)
include_directories (${EXAMPLES_SOURCE_DIR}/../io)
include_directories (${EXAMPLES_SOURCE_DIR}/../gradientdescent)

# Third party internal includes
include_directories (${EXAMPLES_SOURCE_DIR}/../thirdparty)

add_executable (test_fmm test_fmm.cpp
                        ../console/console.cpp
                        ../ndgridmap/cell.cpp
                        ../fmm/fmdata/fmcell.cpp
                        )
target_link_libraries (test_fmm X11 pthread)

#add_executable (test_fm2 test_fm2.cpp
#                        ../console/console.cpp
#                        ../ndgridmap/cell.cpp
#                        ../fmm/fmdata/fmcell.cpp
#                        )
#target_link_libraries (test_fm2 X11 pthread)

add_executable (test_fmm_benchmark test_fmm_benchmark.cpp
                        ../console/console.cpp
                        ../ndgridmap/cell.cpp
                        ../fmm/fmdata/fmcell.cpp
                        )
target_link_libraries (test_fmm_benchmark boost_system boost_filesystem boost_program_options)

add_executable (test_fmm3d test_fmm3d.cpp
                ../console/console.cpp
                ../ndgridmap/cell.cpp
                ../fmm/fmdata/fmcell.cpp
                )
target_link_libraries (test_fmm3d)
add_executable (test_fmm3d_benchmark test_fmm3d_benchmark.cpp
                ../console/console.cpp
                ../ndgridmap/cell.cpp
                ../fmm/fmdata/fmcell.cpp
                )
target_link_libraries (test_fmm3d_benchmark)
