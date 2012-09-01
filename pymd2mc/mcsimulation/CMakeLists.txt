cmake_minimum_required (VERSION 2.6)
option(BUILD_TESTS "Build all of unit tests." FALSE)
option(BUILD_CUDA "Build CUDA version of the simulator." TRUE)
project (mcsim)

#Searching CUDA
if(BUILD_CUDA)
    FIND_PACKAGE(CUDA)
    if(APPLE)
    #make it work on osx
        set (CMAKE_OSX_ARCHITECTURES i386)
        set (CUDA_64_BIT_DEVICE_CODE OFF)
    endif()
endif()


# Set default build type
string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
if(NOT ${CMAKE_BUILD_TYPE})
    set(CMAKE_BUILD_TYPE "debug")
message(STATUS "${PROJECT_NAME}: Build type: ${CMAKE_BUILD_TYPE}")
endif(NOT ${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -fstack-protector-all -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CSS_FLAGS_RELEASE "-O2 -DNO_DEBUG -DNDEBUG -fPIC -rdynamic")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_BINARY_DIR}")
set (BASE_SRC_DIR "${CMAKE_SOURCE_DIR}/src")
if (BUILD_TESTS)
    # setup testing
    enable_testing()
    add_subdirectory(gmock-1.6.0)
    find_package(GTest REQUIRED)
    include_directories(gmock-1.6.0/gtest/include)
    include_directories(gmock-1.6.0/include)
    add_subdirectory(tests)
endif() 

add_subdirectory(src)
