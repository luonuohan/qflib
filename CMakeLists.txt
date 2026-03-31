cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "26.0")
project(qflib)

# This section is only for MSVC, GCC/Clang handle these themselves
# Add these new configurations at the top
if(MSVC)
    # Increase compiler heap limit (300% of default)
    add_compile_options(/Zm300)
    # Allow more sections in object files
    add_compile_options(/bigobj)
    # Optimize memory usage
    add_compile_options(/Zc:inline)
endif()

# project-wide configurations
# sets up where to find Armadillo 
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(THIRDPARTY_DIRECTORY "Address/")
    set(ARMA_VERSION 14.2.2)
elseif(APPLE)
    # On macOS, use Homebrew-installed Armadillo
    set(ARMADILLO_PREFIX "/opt/homebrew/opt/armadillo")
else()
    set(THIRDPARTY_DIRECTORY "$ENV{HOME}/thirdparty")
    set(ARMA_VERSION 14.2.2)
endif()

# set location of artifacts
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# set the debug postfix
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_DEBUG_POSTFIX "-gd")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_DEBUG_POSTFIX "")
endif()

# global compiler and linker option adjustments
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "/EHsc /permissive-")
	# use /Zi for Debug configs
	string(REGEX REPLACE "/Z[iI7]" "" TEMP "${CMAKE_CXX_FLAGS_DEBUG}")
    set(CMAKE_CXX_FLAGS_DEBUG "${TEMP} /Zi")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "-fPIC -Wno-deprecated-declarations -Wno-attributes")
else()
    message(FATAL_ERROR "unknown compiler; only MSVC, GNU, and AppleClang are currently supported" )
endif()

add_subdirectory(qflib)
add_subdirectory(pyqflib)
add_subdirectory(tests)
