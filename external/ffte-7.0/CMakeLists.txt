# Specify the minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Define the project name and specify Fortran as the language
project(ffte LANGUAGES Fortran)

# Set the version number
set(PROJECT_VERSION_MAJOR 7)
set(PROJECT_VERSION_MINOR 0)

# Find all source files with .f extension in the current directory
file(GLOB SOURCE_FILES "*.f")

# Set custom compiler flags based on the compiler
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fomit-frame-pointer")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume dummy_aliases -O3")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Msecond_underscore -Mnoframe")
endif()
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3")

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

# Include directories (if there are any header files located in specific directories)
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
