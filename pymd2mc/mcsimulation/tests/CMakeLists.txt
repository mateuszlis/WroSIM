include_directories(${BASE_SRC_DIR})
add_executable (TriangularLatticeTest testUtils.cpp TriangularLatticeTest.cpp)
add_executable (ProteinTriangularLatticeTest testUtils.cpp ProteinTriangularLatticeTest.cpp)
add_executable (ClustersAnalyzerTest ClustersAnalyzerTest.cpp)
add_executable (MsdEnabledLattExTest MsdEnabledLattExTest.cpp)
#every test has its own valgrind version
#so it is easy to turn off time-consuming valgrind tests in ctest while
#debugging
add_test(TriangularLatticeTest ${EXECUTABLE_OUTPUT_PATH}/TriangularLatticeTest)
add_test(TriangularLatticeTest|Valgrind valgrind --leak-check=full
    --undef-value-errors=no -v --error-exitcode=1 ${EXECUTABLE_OUTPUT_PATH}/TriangularLatticeTest --gtest_filter=-*ConstructException) #exclude exception test, because there is sth wrong with it in OS X
add_test(ProteinTriangularLatticeTest ${EXECUTABLE_OUTPUT_PATH}/ProteinTriangularLatticeTest)
add_test(ProteinTriangularLatticeTest|Valgrind valgrind --leak-check=full
    --undef-value-errors=no -v --error-exitcode=1 ${EXECUTABLE_OUTPUT_PATH}/ProteinTriangularLatticeTest ) #exclude exception test, because there is sth wrong with it in OS X
add_test(ClustersAnalyzerTest  ${EXECUTABLE_OUTPUT_PATH}/ClustersAnalyzerTest)
add_test(ClustersAnalyzerTest|Valgrind valgrind --leak-check=full
    --undef-value-errors=no -v --error-exitcode=1 ${EXECUTABLE_OUTPUT_PATH}/ClustersAnalyzerTest)
add_test(MsdEnabledLattExTest  ${EXECUTABLE_OUTPUT_PATH}/MsdEnabledLattExTest)
add_test(MsdEnabledLattExTest|Valgrind valgrind --leak-check=full
    --undef-value-errors=no -v --error-exitcode=1 ${EXECUTABLE_OUTPUT_PATH}/MsdEnabledLattExTest)
#cppcheck
add_test(cppcheckInfo  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=information -q )
add_test(cppcheckUnusedFunctions  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=unusedFunction -q )
add_test(cppcheckStyle  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=style -q )
add_test(cppcheckAll  cppcheck ${BASE_SRC_DIR} --error-exitcode=1 --enable=all -q )

target_link_libraries(TriangularLatticeTest TriangularLattice gtest gtest_main gmock)
target_link_libraries(ProteinTriangularLatticeTest TriangularLattice ProteinTriangularLattice gtest gtest_main gmock)
target_link_libraries(ClustersAnalyzerTest ClustersAnalyzer gtest gtest_main)
target_link_libraries(MsdEnabledLattExTest TriangularLattice gtest gtest_main)
