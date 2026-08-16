cmake_minimum_required(VERSION 3.15...3.26)
project(fTetWildWrapper_project LANGUAGES CXX)

find_package(Python 3.10
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

# Get nanobind cmake path through python
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
  # name of the extension
  PyfTetWildWrapper

  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  STABLE_ABI

  src/FTetWildWrapper.cpp
)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./src/pytetwild)

# Set the path to the fTetWild project
set(fTetWild_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/fTetWild")

# Enable aggressive optimization. -fsigned-char keeps `char` signed on
# aarch64 (where it defaults to unsigned); fTetWild stores -1 sentinels
# in std::array<char, 4>, which is only valid with signed char.
#
# -ffp-contract=off is geogram's own requirement: exact predicates assume
# every floating point operation rounds, and a fused multiply-add skips
# one, so both its Linux platform files pass this flag and say why. The
# Darwin ones never did, which cost nothing while macOS was x86-64, where
# the baseline ISA has no FMA to contract into, and breaks on arm64, where
# it does. Without it the walk in Delaunay3d::locate_inexact runs to its
# 2500-step cycle guard on nearly every point rather than converging in a
# few, and tetrahedralize runs for minutes instead of seconds. It goes in
# the global flags rather than on geogram because fTetWild carries its own
# copy of Shewchuk's predicates.c, which needs it just as much. See #47.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fsigned-char -ffp-contract=off")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -O3 -fsigned-char -ffp-contract=off")

# Include the fTetWild project as a subdirectory
add_subdirectory(${fTetWild_DIR} EXCLUDE_FROM_ALL)

# Include directories from fTetWild required for the wrapper
target_include_directories(PyfTetWildWrapper PUBLIC ${fTetWild_DIR}/src)

# Link the FloatTetwild library (from fTetWild)
target_link_libraries(PyfTetWildWrapper PRIVATE FloatTetwild)
target_compile_features(PyfTetWildWrapper PUBLIC cxx_std_11)

# still necessary?
if(WIN32)
    set_target_properties(PyfTetWildWrapper PROPERTIES SUFFIX ".pyd")
    foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG_UPPER)
        set_target_properties(PyfTetWildWrapper PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_CURRENT_SOURCE_DIR}/src/pytetwild"
            LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_CURRENT_SOURCE_DIR}/src/pytetwild"
        )
    endforeach()
endif()

install(TARGETS PyfTetWildWrapper LIBRARY DESTINATION pytetwild)
