cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(xyz2pcd)

find_package(PCL 1.2 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable (xyz2pcd xyz2pcd.cpp)
target_link_libraries (xyz2pcd ${PCL_LIBRARIES})

add_executable (xyzrgb2pcd xyzrgb2pcd.cpp)
target_link_libraries (xyzrgb2pcd ${PCL_LIBRARIES})
