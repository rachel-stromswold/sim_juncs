cmake_minimum_required(VERSION 3.3.2)
project(sim_geom)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

#allow the user to specify what type of geometry they'd like to use
if(NOT JUNCTION_TYPE)
    set(JUNCTION_TYPE Au_graphene_box)
endif()

#set names for applications
set(TEST_EXE test)
set(APP_EXE sim_geom)

#find all the libraries we need for building
find_library(meep libmeep.a REQUIRED HINTS $ENV{HOME}/install/lib)
link_directories($ENV{HOME}/install/lib)
#find_library(hdf5 REQUIRED)
#find_library(ctlgeom REQUIRED)

#Cmake can't find the location of eigen on the cluster so we specify that manually
set(Eigen3_DIR "/software/eigen/3.3.7/b1/share/eigen3/cmake")
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

#add targets for libraries
add_executable( ${TEST_EXE} src/cgs.cpp src/disp.cpp src/main_test.cpp )
add_executable( ${APP_EXE} src/cgs.cpp src/disp.cpp src/main.cpp )

#linking
target_link_libraries(${APP_EXE} PRIVATE Eigen3::Eigen meep hdf5 ctlgeom)
target_include_directories( ${APP_EXE} PRIVATE "~/install/include" )
target_include_directories( ${TEST_EXE} PRIVATE "/usr/include/doctest" )
target_link_libraries(${TEST_EXE} PRIVATE Eigen3::Eigen meep)
target_include_directories( ${TEST_EXE} PRIVATE "~/install/include" )
target_include_directories( ${TEST_EXE} PRIVATE "/usr/include/doctest" )

#copy files needed for unit testing
configure_file(tests/test.eps test.eps COPYONLY)
configure_file(tests/test.conf test.conf COPYONLY)

#copy the junction geometry files
configure_file(junctions/${JUNCTION_TYPE}/junc_template.geom COPYONLY)
configure_file(junctions/${JUNCTION_TYPE}/params.conf COPYONLY)

#copy utility scripts that are used
configure_file(scripts/run.sh COPYONLY)
configure_file(scripts/utils.py COPYONLY)
configure_file(scripts/set_dims.py COPYONLY)
configure_file(scripts/check_enes.py COPYONLY)
