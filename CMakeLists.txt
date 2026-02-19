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

# Enable aggressive optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -O3")

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
