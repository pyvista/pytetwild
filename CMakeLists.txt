cmake_minimum_required(VERSION 3.15...3.26)
project(fTetWildWrapper)

# Set the path to the fTetWild project
set(fTetWild_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/fTetWild")

# Add the pybind11 submodule
set(PYBIND11_NEWPYTHON ON)
set(pybind11_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/pybind11")
add_subdirectory(${pybind11_DIR})

# Include the fTetWild project as a subdirectory
add_subdirectory(${fTetWild_DIR} EXCLUDE_FROM_ALL)

pybind11_add_module(PyfTetWildWrapper MODULE "${CMAKE_CURRENT_SOURCE_DIR}/src/FTetWildWrapper.cpp")

# Include directories from fTetWild required for the wrapper
target_include_directories(PyfTetWildWrapper PUBLIC ${fTetWild_DIR}/src)

# Link the FloatTetwild library (from fTetWild)
target_link_libraries(PyfTetWildWrapper PRIVATE FloatTetwild)
target_compile_features(PyfTetWildWrapper PUBLIC cxx_std_11)

install(TARGETS PyfTetWildWrapper LIBRARY DESTINATION ./pytetwild)
