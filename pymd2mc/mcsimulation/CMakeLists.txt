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

add_library(TriangularLattice TriangularLattice.cpp)
add_executable(mcsim mcsim.cpp)
target_link_libraries(mcsim TriangularLattice)