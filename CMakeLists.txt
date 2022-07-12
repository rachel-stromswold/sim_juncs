cmake_minimum_required(VERSION 3.3.2)
project(sim_geom)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

set(TEST_EXE test)
set(APP_EXE sim_geom)

#find_library(meep REQUIRED)
#find_library(hdf5 REQUIRED)
#find_library(ctlgeom REQUIRED)

set(Eigen3_DIR "/software/eigen/3.3.7/b1/share/eigen3/cmake")
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable( ${TEST_EXE} src/cgs.cpp src/disp.cpp src/main_test.cpp )
add_executable( ${APP_EXE} src/cgs.cpp src/disp.cpp src/main.cpp )

target_link_libraries(${TEST_EXE} PUBLIC Eigen3::Eigen meep)
target_include_directories( ${TEST_EXE} PRIVATE "~/local_include/include" )
target_include_directories( ${TEST_EXE} PRIVATE "/usr/include/doctest" )
target_link_libraries(${APP_EXE} PUBLIC Eigen3::Eigen meep hdf5 ctlgeom)

configure_file(tests/test.mt test.mt COPYONLY)
configure_file(tests/params.conf params.conf COPYONLY)
