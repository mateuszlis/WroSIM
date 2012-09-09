project(cuda_metropolis)
message(STATUS "==============================CUDA enabled! ==================================")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/;${CMAKE_MODULE_PATH}")
message(STATUS "Module path is now set to: " ${CMAKE_MODULE_PATH} )

set(COMMON_SOURCE_DIR "${PROJECT_SOURCE_DIR}")
message(STATUS "CommonSource Dir is set to: ${COMMON_SOURCE_DIR} ")

#Searching CUDA
FIND_PACKAGE(CUDA)

# Set the flags
SET(CUDA_PROPAGATE_HOST_FLAGS ON CACHE BOOL "Propagate C/CXX Flags and friends to the host compiler in NVCC via -Xompile  " FORCE)

include(PrintListMacro)
include(FindCUDASDK)


#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#Append includes
string(REPLACE "/C" "/" NVIDIA_SDK_DIR "${CUDA_SDK_ROOT_DIR}")
message(STATUS "NVIDIA Sample ROOT Dir is: ${NVIDIA_SDK_DIR}")

#include directories
set(INCLUDE_DIRS 
"${CUDA_CUT_INCLUDE_DIR}"
"${NVIDIA_SDK_DIR}/shared/inc"
"${NVIDIA_SDK_DIR}CUDALibraries/common/UtilNPP"
"${NVIDIA_SDK_DIR}CUDALibraries/common/FreeImage/include"
"${PROJECT_SOURCE_DIR}/inc"
"${COMMON_SOURCE_DIR}/inc"
)
include_directories(${INCLUDE_DIRS})
PRINTLIST("Include Directories are:" "${INCLUDE_DIRS}")

set( CUDA_NPP_LIBRARY
    /usr/local/cuda/lib/libnpp.dylib )
set( CUDA_NPP_UTIL_LIBRARY
    ${NVIDIA_SDK_DIR}CUDALibraries/common/lib/libUtilNPP_i386.a )
set( CUDA_FREEIMAGE_LIBRARY
    ${NVIDIA_SDK_DIR}CUDALibraries/common/FreeImage/lib/darwin/libfreeimage_10_6.a
    )
set(LINKED_LIBS  ${CUDA_CUDA_LIBRARY} ${CUDA_CUT_LIBRARY} ${CUDA_NPP_LIBRARY} ${CUDA_FREEIMAGE_LIBRARY}
    ${CUDA_NPP_UTIL_LIBRARY} TriangularLattice Metropolis )
PRINTLIST("Linked libraries are:" "${LINKED_LIBS}")

SET(CMAKE_DEBUG_POSTFIX "D")

set_source_files_properties( INCLUDE_FILES PROPERTIES HEADER_FILE_ONLY TRUE)


message( STATUS ${COMMON_FILES_SRC})
#The rule to build the executable program
CUDA_ADD_LIBRARY(Metropolis_cuda  Metropolis.cu )
target_link_libraries(Metropolis_cuda  ${LINKED_LIBS} )
message(STATUS "=========================================================================")
