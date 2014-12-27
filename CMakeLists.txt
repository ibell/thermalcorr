cmake_minimum_required(VERSION 2.8)

# Project name
set(project_name "ThermalCorr")
set(app_name ${project_name})
project(${project_name})

file(GLOB_RECURSE APP_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB COOLPROP_SOURCES "${CMAKE_SOURCE_DIR}/externals/CoolProp/CoolProp/*.cpp")

list (APPEND APP_SOURCES ${COOLPROP_SOURCES})
list (APPEND APP_SOURCES "${CMAKE_SOURCE_DIR}/main.cpp")

include_directories("${CMAKE_SOURCE_DIR}/externals/coolprop" )

add_executable(ThermalCorr ${APP_SOURCES})