project(cuda_metropolis)
message(STATUS "==============================CUDA enabled! ==================================")

set (CMAKE_OSX_ARCHITECTURES i386)
set (CUDA_64_BIT_DEVICE_CODE OFF)
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
#"${NVIDIA_SDK_DIR}/C/common/inc"
"${CUDA_CUT_INCLUDE_DIR}"
"${NVIDIA_SDK_DIR}/shared/inc"
"${PROJECT_SOURCE_DIR}/inc"
"${COMMON_SOURCE_DIR}/inc"
)
include_directories(${INCLUDE_DIRS})
PRINTLIST("Include Directories are:" "${INCLUDE_DIRS}")

##link directories
#set(LINK_DIRS
##"${NVIDIA_SDK_DIR}/C/common/lib"
#"${CUDA_CUT_LIBRARY}"
#"${NVIDIA_SDK_DIR}/shared/lib"
#)
##link_directories(${LINK_DIRS})

#PRINTLIST("Link directories are:" "${LINK_DIRS}")
set(LINKED_LIBS  ${CUDA_CUDA_LIBRARY} ${CUDA_CUT_LIBRARY} TriangularLattice)
PRINTLIST("Linked libraries are:" "${LINKED_LIBS}")



set(SOURCE_FILES 

)

set(INCLUDE_FILES

) 

SET(CMAKE_DEBUG_POSTFIX "D")

set_source_files_properties( INCLUDE_FILES PROPERTIES HEADER_FILE_ONLY TRUE)


message( STATUS ${COMMON_FILES_SRC})
#The rule to build the executable program
CUDA_ADD_LIBRARY(Metropolis_cuda  Metropolis.cu )
target_link_libraries(Metropolis_cuda  ${LINKED_LIBS})

message(STATUS "=========================================================================")
