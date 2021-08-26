# stop if cmake version below 3.10
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# project name
project(ddcMD)

# project version
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(VERSION_STRING ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# enable fortran, c, and c++ language
enable_language(C CXX)
message(STATUS "Compiler is ${CMAKE_CXX_COMPILER_ID}")

option(USE_GPU  "Build on GPU" ON)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

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
    message(FATAL_ERROR "ddcMD does not support a ${CMAKE_SYSTEM_NAME} system")
endif()

cmake_host_system_information(RESULT host_name_raw QUERY HOSTNAME)
string(REGEX MATCH "^([^.]+)" host_name ${host_name_raw})
string(REGEX REPLACE "[0-9]" "" host_name "${host_name}")
message(STATUS "Machine name is " ${host_name})

#if (CMAKE_CXX_COMPILER_ID STREQUAL "XL")
#    if(${CMAKE_VERSION} VERSION_GREATER "3.17.999")
#    message(FATAL_ERROR "For XL compiler please use CMAKE version less than 3.17")
#    endif()
#endif()

# find OpenMP library
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    add_definitions(-DWITH_OMP)
    include_directories( ${OpenMP_C_INCLUDE_DIRS} )
endif(OpenMP_FOUND)

# find MPI library
find_package(MPI REQUIRED COMPONENTS C CXX MPICXX)
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

# Use emplicit list of source files
file(GLOB cpu_source_files "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB gpu_source_files "${CMAKE_SOURCE_DIR}/src/*.cu")

# Remove explicit list of source file treee sourceFileTree.cmake
# defined in ${CMAKE_CURRENT_SOURCE_DIR}/cmake
# include(sourceFileTree)

if(USE_GPU)
    add_definitions(-DUSE_GPU=1)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_RUNTIME_LIBRARY Shared)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    message(STATUS "CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES=${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    include_directories(${CMAKE_SOURCE_DIR}/cub)
    list(REMOVE_ITEM cpu_source_files "${CMAKE_SOURCE_DIR}/src/constructList.c")
    add_executable(ddcMD-gpu ${cpu_source_files} ${gpu_source_files})
    set_target_properties(ddcMD-gpu PROPERTIES OUTPUT_NAME "ddcMD_GPU_${host_name}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_BUILD_TYPE}")
    #Note that only support P100 and V100, have not tested against other GPUs
    if(${CMAKE_VERSION} VERSION_GREATER "3.17.999")
        set_property(TARGET ddcMD-gpu  PROPERTY CUDA_ARCHITECTURES 70 60)
    else()
        target_compile_options(ddcMD-gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_70,code=sm_70 -gencode=arch=compute_60,code=sm_60>)
    endif()
    target_compile_options(ddcMD-gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-g>)
    #link_directories("/usr/tce/packages/cuda/cuda-10.1.243/lib64")
    #target_compile_options(ddcMD-gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_70,code=sm_70>)
    #target_link_libraries(ddcMD-gpu ${MPI_C_LIBRARIES} ${FFTW_LIBRARIES} ${CUDA_LIBRARIES} )
    #set_target_properties(ddcMD-gpu PROPERTIES CUDA_USE_STATIC_CUDA_RUNTIME OFF)
    #target_link_libraries(ddcMD-gpu ${MPI_C_LIBRARIES} ${FFTW_LIBRARIES} cudart cuda curand nvToolsExt nvrtc)
    target_link_libraries(ddcMD-gpu ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OpenMP_C_LIBRARIES})
    target_link_libraries(ddcMD-gpu fftw3 cudart cuda curand nvToolsExt nvrtc)
    set_target_properties(ddcMD-gpu PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS ddcMD-gpu DESTINATION bin)
else(USE_GPU)
    add_executable(ddcMD-cpu ${cpu_source_files})
    set_target_properties(ddcMD-cpu PROPERTIES OUTPUT_NAME "ddcMD_CPU_${host_name}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_BUILD_TYPE}")
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        target_link_libraries(ddcMD-cpu ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OpenMP_C_LIBRARIES} ${FFTW_C_LIBRARIES})
    else()
        target_link_libraries(ddcMD-cpu ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OpenMP_C_LIBRARIES} fftw3)
    endif()
    set_target_properties(ddcMD-cpu PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
    install(TARGETS ddcMD-cpu DESTINATION bin)
endif(USE_GPU)

