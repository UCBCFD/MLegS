# Specify the minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Define the project name and specify Fortran as the language
project(fm LANGUAGES Fortran)

# Set the version number
set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 4)

# Set the module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/inc)

# Find all source files with .f90 extension in the current directory
file(GLOB SOURCE_FILES "*.f90")

# Set the library name
set(LIB_NAME ${PROJECT_NAME})

# Create a static library from the source files
add_library(${LIB_NAME} STATIC ${SOURCE_FILES})

# Specify the output directory for static libraries
set_target_properties(${LIB_NAME} PROPERTIES
	ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Specify the library version
set_target_properties(${LIB_NAME} PROPERTIES VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

# Include directories/libraries
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_Fortran_MODULE_DIR})

# Custom target to clean up .mod files
add_custom_target(clean_mod_files
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_Fortran_MODULE_DIR}/*.mod
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_Fortran_MODULE_DIR}/*.mod
    COMMENT "Removing all .mod files"
)
