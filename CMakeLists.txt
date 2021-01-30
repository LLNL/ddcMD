# stop if cmake version below 3.8
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# project name
project(ddcMD)

# project version
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(VERSION_STRING ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# enable fortran, c, and c++ language
enable_language(C CXX)

option(USE_GPU  "Build on GPU" ON)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

#set (CMAKE_CXX_STANDARD 11)

# location of additional cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# detect operating system
message(STATUS "Building on a ${CMAKE_SYSTEM_NAME} system")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYSTEM_LINUX -DWITH_MPI -DWITH_PIO -D_GNU_SOURCE)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYSTEM_DARWIN -DWITH_MPI -DWITH_PIO -D_GNU_SOURCE)
    #add_compile_options(-Xpreprocessor -fopenmp)
else()
    message("ddcMD does not support a ${CMAKE_SYSTEM_NAME} system")
endif()

# find MPI library
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# find FFTW library
find_package(FFTW REQUIRED)
include_directories( ${FFTW_INCLUDE_DIRS} )

# find GSL library
find_package(GSL QUIET)
if(GSL_FOUND)
    include_directories( ${GSL_INCLUDE_DIRS} )
    add_definitions(-DHAVE_GSL)
endif(GSL_FOUND)

file(GLOB cpu_source_files "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB gpu_source_files "${CMAKE_SOURCE_DIR}/src/*.cu")

if(USE_GPU)
    add_definitions(-DUSE_GPU=1)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    include_directories(${CMAKE_SOURCE_DIR}/cub)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    list(REMOVE_ITEM cpu_source_files "${CMAKE_SOURCE_DIR}/src/constructList.c")
    add_executable(ddcMD-gpu ${cpu_source_files} ${gpu_source_files})
    set_property(TARGET ddcMD-gpu  PROPERTY CUDA_ARCHITECTURES 70)
    #target_compile_options(ddcMD-gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-std=c++11>) 
    #target_compile_options(ddcMD-gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_70,code=sm_70>)
    target_link_libraries(ddcMD-gpu ${MPI_C_LIBRARIES} ${FFTW_LIBRARIES} cudart cuda curand nvToolsExt nvrtc)
    set_target_properties(ddcMD-gpu PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS ddcMD-gpu DESTINATION bin)
else(USE_GPU)
    add_executable(ddcMD-cpu ${cpu_source_files})
    target_link_libraries(ddcMD-cpu ${MPI_C_LIBRARIES} ${FFTW_C_LIBRARIES})
    set_target_properties(ddcMD-cpu PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS ddcMD-cpu DESTINATION bin)
endif(USE_GPU)

