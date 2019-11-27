cmake_minimum_required(VERSION 3.14)
project (MAT VERSION 0.0.1)

add_definitions(-DGF_PLUGIN_NAME=\"${PROJECT_NAME}\")

find_package(geoflow REQUIRED)

# dependencies
find_package(glm)
find_package(CGAL QUIET COMPONENTS Core)

if(EXISTS "${PROJECT_SOURCE_DIR}/.gitmodules")
  execute_process(
    COMMAND git submodule update --init masbcpp region-grower
    WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/external
    )
endif()

add_subdirectory(external/masbcpp)
include_directories(
  external/masbcpp/src
  external/masbcpp/thirdparty
  external/masbcpp/thirdparty/vrui
  external/region-grower
)

add_library(gfp_masb MODULE
  plugin.cpp
  masb_nodes.cpp
)
target_link_libraries(gfp_masb 
  geoflow-core 
  glm 
  masbcpp
CGAL::CGAL CGAL::CGAL_Core
)

set_target_properties(
  gfp_masb
  PROPERTIES CXX_STANDARD 17
)

install(TARGETS gfp_masb DESTINATION ${geoflow_PLUGIN_FOLDER})