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

# Performance optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG -ffast-math")
    
    # Enable SIMD optimizations based on architecture
    if(APPLE)
        # For universal binaries, use runtime detection instead of compile-time flags
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")  # Let runtime code handle SIMD detection
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=x86-64-v2 -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
    endif()
elseif(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG /fp:fast")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2 /Zi /DNDEBUG /fp:fast")
    
    # Enable SIMD for x86_64 on MSVC
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
    endif()
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
