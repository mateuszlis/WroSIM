include_directories(${BASE_SRC_DIR})
add_executable (TriangularLatticeTest TriangularLatticeTest.cpp)
add_executable (ClustersAnalyzerTest ClustersAnalyzerTest.cpp)

add_test(TriangularLatticeTest ${EXECUTABLE_OUTPUT_PATH}/TriangularLatticeTest)
add_test(ClustersAnalyzerTest  ${EXECUTABLE_OUTPUT_PATH}/ClustersAnalyzerTest)
#cppcheck
add_test(cppcheckInfo  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=information -q )
add_test(cppcheckUnusedFunctions  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=unusedFunction -q )
add_test(cppcheckStyle  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=style -q )

target_link_libraries(TriangularLatticeTest TriangularLattice)
target_link_libraries(ClustersAnalyzerTest ClustersAnalyzer)
