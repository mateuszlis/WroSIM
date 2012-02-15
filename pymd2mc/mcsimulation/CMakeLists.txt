cmake_minimum_required (VERSION 2.6)
option(BUILD_TESTS "Build all of foo's unit tests." FALSE)
project (mcsim)

#Searching CUDA
FIND_PACKAGE(CUDA)

#make it work on osx
set (CMAKE_OSX_ARCHITECTURES i386)
set (CUDA_64_BIT_DEVICE_CODE OFF)

# Set default build type
string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
message(STATUS "${PROJECT_NAME}: Build type: ${CMAKE_BUILD_TYPE}")
if(NOT ${CMAKE_BUILD_TYPE})
    set(CMAKE_BUILD_TYPE "debug")
endif(NOT ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -fstack-protector-all -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CSS_FLAGS_RELEASE "-O2 -DNO_DEBUG -DNDEBUG -fPIC -rdynamic")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_BINARY_DIR}")
set (BASE_SRC_DIR ${CMAKE_SOURCE_DIR})
if (BUILD_TESTS)
    # setup testing
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIR})
    include_directories(${GMOCK_INCLUDE_DIR})
    add_subdirectory(tests)
endif() 


add_library(mcsimCommandLine mcsimCommandLine.cxx)
add_library(analyzePyrkovaCommandLine analyzePyrkovaCommandLine.cxx)
add_library(TriangularLattice TriangularLattice.cpp)

add_subdirectory(gpu)

add_library(Metropolis Metropolis.cpp)
add_library(ClustersAnalyzer ClustersAnalyzer.cpp)
add_library(Atom Atom.cpp)

add_executable(mcsim mcsim.cpp)
add_executable(analyze analyze.cpp)
add_executable(analyzePyrkova analyzePyrkova.cpp)
if(CUDA_FOUND)
    message(STATUS "CUDA")
    target_link_libraries(mcsim TriangularLattice mcsimCommandLine Metropolis_cuda)
else(CUDA_FOUND)
    target_link_libraries(mcsim TriangularLattice mcsimCommandLine Metropolis)
endif(CUDA_FOUND)
target_link_libraries(analyzePyrkova ClustersAnalyzer Atom
    analyzePyrkovaCommandLine)

