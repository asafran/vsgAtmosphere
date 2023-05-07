cmake_minimum_required(VERSION 3.7)

project(vsgAtmosphere
    VERSION 0.0.1
    DESCRIPTION "Brunton's sky model port to VSG"
    LANGUAGES CXX C
)

# build all examples into the bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(VSGATMOSPHERE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Root source directory of vsgXchange.")
set(VSGATMOSPHERE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE INTERNAL "Root binary directory of vsgXchange.")

# set the use of C++17 globally as all examples require it
set(CMAKE_CXX_STANDARD 17)

# Find Vulkan and the VSG
if (VULKAN_SDK)
    set(ENV{VULKAN_SDK} ${VULKAN_SDK})
endif()

set(VSG_MIN_VERSION 1.0.3)
find_package(vsg ${VSG_MIN_VERSION})

vsg_setup_dir_vars()
vsg_setup_build_vars()

# find the optional vsgExamples that can be used for reading and range of image and 3d model formats and shader compilation
find_package(vsgXchange 1.0.2 QUIET)

# for generated cmake support files
set(FIND_DEPENDENCY ${FIND_DEPENDENCY} "find_dependency(vsg ${VSG_MIN_VERSION} REQUIRED)")

#vsg_add_target_uninstall()

# source directory for main vsgXchange library
add_subdirectory(src)
add_subdirectory(examples/vsgatmosphere)

vsg_add_feature_summary()
