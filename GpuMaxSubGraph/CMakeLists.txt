cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(subgraph_search VERSION 5.0 LANGUAGES CXX CUDA)

# ########################################################################
# ### Tweak the default configuration
# ########################################################################

# define an option for using native cpp code
set(TARGET_CPU_ARCHITECTURE "generic" CACHE STRING "The target cpu for architecture-specific optimizations")
set_property(CACHE TARGET_CPU_ARCHITECTURE PROPERTY STRINGS generic native)


# the idea is that we want to build the optimized version, unless the user says otherwise
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "The type of build" FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()


# export also compile commands to have a better integration with the ide
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



# enable colored compiler output (with Ninja)
if("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-fdiagnostics-color=always)
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-fcolor-diagnostics)
  endif()
endif()

# ########################################################################
# ### Add external dependency
# ########################################################################
find_program(VALGRIND_EXECUTABLE valgrind)

if(VALGRIND_EXECUTABLE)
  set(MEMORYCHECK_COMMAND "${VALGRIND_EXECUTABLE}" "--leak-check=full" "--show-leak-kinds=all" "--track-origins=yes" "--verbose" "--log-file=valgrind.log")
  set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind.log")
  add_custom_target(memcheck COMMAND ${MEMORYCHECK_COMMAND} $<TARGET_FILE: main>)
endif()

# ########################################################################
# ### Define the application sources
# ########################################################################

# application headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND header_files
  "${header_path}/main.hpp"
  "${header_path}/cuda_header.h"
  "${header_path}/gpu_header.hpp"
)

# application sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files
  "${source_path}/Label.cpp"
  "${source_path}/gen_rotations.cpp"
)

# ########################################################################
# ### Define the compilation step
# ########################################################################

# define the target and its properties

add_executable(main ${header_files} "${source_path}/main.cu"
        "${source_path}/gen_bond_labels.cpp"
        "${source_path}/gen_rotations.cpp"
        "${source_path}/select_vertex.cpp"
        "${source_path}/hood.cpp"
        "${source_path}/Label.cpp"
        "${source_path}/smiles_mcs.cpp"
        "${source_path}/calc_bound.cpp"
        "${source_path}/select_label.cpp"
        "${source_path}/gen_initial_labels.cpp"
        "${source_path}/gen_ring_classes.cpp"
        "${source_path}/gen_rotations.cpp"
        "${source_path}/search_mcs.cpp"
        "${source_path}/mol_mcs.cpp"
        "${source_path}/mc_split.cpp"
        "${source_path}/g2mol.cpp"
        "${source_path}/gpu_mc_split.cpp"
        "${source_path}/pair_vertex.cpp"
        "${source_path}/kernel.cu"
        "${source_path}/initialized.cpp"
        "${source_path}/mcs_iterative.cpp"
)
set_target_properties(main PROPERTIES
    COMPILE_DEFINITIONS "TESTING"
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)



# enable link-time optimizations
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)



# since RDKit doesn't handle the dependencies in a correct way, we need to improvise
# NOTE: by using this way of importing RDKit, we need to manually include its dependencies
# Get the installation prefix of RDKit

get_filename_component(RDKIT_INSTALL_PREFIX "/home/francesco/miniconda3/envs/my-rdkit-env" REALPATH)

# Add the installation prefix of RDKit to CMAKE_PREFIX_PATH

list(APPEND CMAKE_PREFIX_PATH "${RDKIT_INSTALL_PREFIX}")

# the rdkit to perform the heavy lifting
find_package(RDKit REQUIRED)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

get_target_property(RDKIT_LIB_FULLPATH RDKit::RDGeneral LOCATION)
get_filename_component(RDKIT_LIB_DIRPATH "${RDKIT_LIB_FULLPATH}" DIRECTORY)
cmake_path(GET RDKIT_LIB_DIRPATH PARENT_PATH RDKIT_INSTALL_PREFIX)
include_directories(${RDKit_INCLUDE_DIRS}) #modified

include_directories(main PUBLIC "${RDKIT_INSTALL_PREFIX}/include/rdkit")
target_link_directories(main PUBLIC "${RDKIT_INSTALL_PREFIX}/lib")
target_link_directories(main PUBLIC "${RDKIT_INSTALL_PREFIX}/lib64")
target_link_libraries(main PUBLIC
        -lRDKitFileParsers
        -lRDKitGraphMol
        -lRDKitRDGeneral
        -lRDKitSmilesParse
        -lRDKitSubstructMatch
)



message(STATUS "RDKit Include Directories: ${RDKit_INCLUDE_DIRS}")
message(STATUS "RDKit Library Directories: ${RDKit_LIBRARY_DIRS}")

# CUDA specific settings
find_package(CUDA REQUIRED)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Set CUDA architectures to compile for
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 60 61 70 75 80)
endif()


