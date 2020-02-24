cmake_minimum_required(VERSION 3.14)
project (MAT VERSION 0.0.1)

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

set(GF_PLUGIN_NAME ${PROJECT_NAME})
set(GF_PLUGIN_TARGET_NAME "gfp_masb")
set(GF_PLUGIN_REGISTER ${PROJECT_SOURCE_DIR}/register.hpp)
geoflow_create_plugin(masb_nodes.cpp)

target_link_libraries( gfp_masb 
  geoflow-core 
  glm 
  masbcpp
  CGAL::CGAL CGAL::CGAL_Core
)

if(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
  target_compile_options(gfp_masb PRIVATE "/Zc:__cplusplus")
endif()
