cmake_minimum_required(VERSION 3.15)

project(evans_comb_filter
    VERSION 1.1.2
    LANGUAGES C CXX
)

# ---- language standards ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---- macOS universal binary support ----
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "" FORCE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "" FORCE)
endif()

# ---- sane default build type ----
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

# ---- compiler flags (DSP-safe) ----
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -O2
        -g
        -fno-fast-math
    )
elseif(MSVC)
    add_compile_options(/O2 /Zi)
endif()

# ---- fetch CLAP headers ----
include(FetchContent)

FetchContent_Declare(
    clap
    GIT_REPOSITORY https://github.com/free-audio/clap.git
    GIT_TAG        1.2.6
)

FetchContent_MakeAvailable(clap)

# ---- plugin target ----
add_library(evans_comb_filter SHARED
    plugin.h
    plugin.cpp
    dsp.cpp
    filters.cpp
    params.cpp
    utils.cpp
)

# ---- include paths ----
target_include_directories(evans_comb_filter PRIVATE
    ${clap_SOURCE_DIR}/include
)

# ---- CLAP plugin output settings ----
set_target_properties(evans_comb_filter PROPERTIES
    PREFIX ""
    SUFFIX ".clap"
    OUTPUT_NAME "EvansHarmonicIsolator"
)

# ---- platform sanity ----
if(APPLE)
    set_target_properties(evans_comb_filter PROPERTIES
        MACOSX_RPATH ON
    )
endif()
