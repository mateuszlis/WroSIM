cmake_minimum_required (VERSION 2.6)
project (mcsim)
enable_testing()
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
set (directories
    tests
    )
foreach (f ${directories})
    add_subdirectory(${f})
endforeach(f)
add_library(mcsimCommandLine mcsimCommandLine.cxx)
add_library(analyzePyrkovaCommandLine analyzePyrkovaCommandLine.cxx)
add_library(TriangularLattice TriangularLattice.cpp)
add_library(Metropolis Metropolis.cpp)
add_library(ClustersAnalyzer ClustersAnalyzer.cpp)
add_library(Atom Atom.cpp)

add_executable(mcsim mcsim.cpp)
add_executable(analyze analyze.cpp)
add_executable(analyzePyrkova analyzePyrkova.cpp)

target_link_libraries(mcsim TriangularLattice mcsimCommandLine Metropolis)
target_link_libraries(analyzePyrkova ClustersAnalyzer Atom
    analyzePyrkovaCommandLine)
