cmake_minimum_required(VERSION 3.15...3.26)
project(fTetWildWrapper_project LANGUAGES CXX)

find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./src/pytetwild)

if($ENV{USE_MAVX})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
endif()

# Set the path to the fTetWild project
set(fTetWild_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/fTetWild")

# Enable aggressive optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -O3")

# # Add the pybind11 submodule
# set(PYBIND11_NEWPYTHON ON)
# set(pybind11_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/pybind11")
# add_subdirectory(${pybind11_DIR})


# Include the fTetWild project as a subdirectory
add_subdirectory(${fTetWild_DIR} EXCLUDE_FROM_ALL)

# pybind11_add_module(PyfTetWildWrapper MODULE "${CMAKE_CURRENT_SOURCE_DIR}/src/FTetWildWrapper.cpp")

nanobind_add_module(
  # name of the extension
  PyfTetWildWrapper

  # Target the stable ABI for Python 3.12+, which reduces
  # the number of binary wheels that must be built. This
  # does nothing on older Python versions
  STABLE_ABI

  "${CMAKE_CURRENT_SOURCE_DIR}/src/FTetWildWrapper.cpp"
)


# Include directories from fTetWild required for the wrapper
target_include_directories(PyfTetWildWrapper PUBLIC ${fTetWild_DIR}/src)

# Link the FloatTetwild library (from fTetWild)
target_link_libraries(PyfTetWildWrapper PRIVATE FloatTetwild)
target_compile_features(PyfTetWildWrapper PUBLIC cxx_std_11)

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

install(TARGETS PyfTetWildWrapper LIBRARY DESTINATION ./pytetwild)
