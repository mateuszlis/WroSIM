include_directories(${BASE_SRC_DIR})
add_executable (TriangularLatticeTest TriangularLatticeTest.cpp)
add_test(TriangularLatticeTest ${EXECUTABLE_OUTPUT_PATH}/TriangularLatticeTest)
target_link_libraries(TriangularLatticeTest TriangularLattice)