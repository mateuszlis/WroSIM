cmake_minimum_required (VERSION 2.6)
project (mcsim)
enable_testing()
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
