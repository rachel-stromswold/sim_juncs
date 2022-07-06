cmake_minimum_required(VERSION 2.8.9)
project(sim_geom)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

set(TEST_EXE test)
set(APP_EXE sim_geom)

include_directories( "src" "/usr/include/doctest" )
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
#find_library(meep REQUIRED)
#find_library(hdf5 REQUIRED)
#find_library(ctlgeom REQUIRED)

add_executable( ${TEST_EXE} src/cgs.cpp src/main_test.cpp )
add_executable( ${APP_EXE} src/cgs.cpp src/disp.cpp src/main.cpp )

target_link_libraries(${TEST_EXE} PUBLIC Eigen3::Eigen)
target_link_libraries(${APP_EXE} PUBLIC Eigen3::Eigen meep hdf5 ctlgeom)

configure_file(tests/test.mt test.mt COPYONLY)
