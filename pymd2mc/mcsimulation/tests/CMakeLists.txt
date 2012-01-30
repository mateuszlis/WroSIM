include_directories(${BASE_SRC_DIR})
add_executable (TriangularLatticeTest TriangularLatticeTest.cpp)
add_executable (ClustersAnalyzerTest ClustersAnalyzerTest.cpp)
#every test has its own valgrind version
#so it is easy to turn off time-consuming valgrind tests in ctest while
#debugging
add_test(TriangularLatticeTest ${EXECUTABLE_OUTPUT_PATH}/TriangularLatticeTest)
add_test(TriangularLatticeTest|Valgrind valgrind --leak-check=full
    --undef-value-errors=no -v --error-exitcode=1 ${EXECUTABLE_OUTPUT_PATH}/TriangularLatticeTest)
add_test(ClustersAnalyzerTest  ${EXECUTABLE_OUTPUT_PATH}/ClustersAnalyzerTest)
add_test(ClustersAnalyzerTest|Valgrind valgrind --leak-check=full
    --undef-value-errors=no -v --error-exitcode=1 ${EXECUTABLE_OUTPUT_PATH}/ClustersAnalyzerTest)
#cppcheck
add_test(cppcheckInfo  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=information -q )
add_test(cppcheckUnusedFunctions  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=unusedFunction -q )
add_test(cppcheckStyle  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=style -q )

target_link_libraries(TriangularLatticeTest TriangularLattice gtest gtest_main
    gmock)
target_link_libraries(ClustersAnalyzerTest ClustersAnalyzer gtest gtest_main)
