cmake_minimum_required(VERSION 2.8)

###     FLUIDS, MIXTURES JSON       ###
find_package (PythonInterp 2.7 REQUIRED)
add_custom_target(generate_headers
                  COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_SOURCE_DIR}/externals/coolprop/dev/generate_headers.py")
                  
# Project name
set(project_name "ThermalCorr")
set(app_name ${project_name})
project(${project_name})

file(GLOB THERMALCORR_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE COOLPROP_SOURCES "${CMAKE_SOURCE_DIR}/externals/coolprop/src/*.cpp")

include_directories("${CMAKE_SOURCE_DIR}/externals/coolprop/include" "${CMAKE_SOURCE_DIR}/externals/coolprop/externals/Eigen" "${CMAKE_SOURCE_DIR}/externals/coolprop/externals/REFPROP-headers")

add_library(ThermalCorr STATIC ${THERMALCORR_SOURCES})
add_library(CoolProp STATIC ${COOLPROP_SOURCES})

add_dependencies (CoolProp generate_headers)

add_executable(main "${CMAKE_SOURCE_DIR}/main.cpp")
target_link_libraries(main "$<TARGET_FILE:ThermalCorr>" "$<TARGET_FILE:CoolProp>")