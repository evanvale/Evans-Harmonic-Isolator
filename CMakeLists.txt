cmake_minimum_required(VERSION 3.15)
project(evans_comb_filter)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Apple universal binary support
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64")  # Universal binary for both Intel and Apple Silicon
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")     # Support macOS Catalina and newer
endif()

# Try gentle optimization that preserves audio quality
# Use -O1 (basic optimization) without aggressive math changes
set(CMAKE_CXX_FLAGS "-g -O1")
set(CMAKE_C_FLAGS "-g -O1")

# Override any build type flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-g -O1 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O1 -DNDEBUG")
    
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O1 /Zi /DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O1 /Zi /DNDEBUG")
endif()

# Fetch CLAP headers from official repository
include(FetchContent)

FetchContent_Declare(
  clap
  GIT_REPOSITORY https://github.com/free-audio/clap.git
  GIT_TAG        1.2.6
)

FetchContent_MakeAvailable(clap)

# create the plugin library
add_library(evans_comb_filter SHARED
    plugin.h
    plugin.cpp
    dsp.cpp
    filters.cpp
    params.cpp
    utils.cpp
)

# Link CLAP headers
target_include_directories(evans_comb_filter PRIVATE ${clap_SOURCE_DIR}/include)

# set proper library properties for CLAP plugin
set_target_properties(evans_comb_filter PROPERTIES
    PREFIX ""
    SUFFIX ".clap"
    OUTPUT_NAME "EvansHarmonicIsolator"
)
