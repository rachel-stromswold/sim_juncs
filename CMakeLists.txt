cmake_minimum_required(VERSION 3.3.2)
project(sim_geom)
set (CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

#set names for applications
if(CMAKE_BUILD_TYPE MATCHES "^[Dd]ebug")
    set(TEST_EXE test)
endif()
set(APP_EXE sim_geom)

#find all the libraries we need for building
find_library(meep libmeep.a REQUIRED HINTS $ENV{HOME}/install/lib)
#we need to install cmake to a custom location
#set(HDF5_DIR "/home/sstromsw/install/cmake")
find_package(HDF5 REQUIRED COMPONENTS CXX HL)

find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)

#Cmake can't find the location of eigen on the cluster so we specify that manually
set(Eigen3_DIR "/software/eigen/3.3.7/b1/share/eigen3/cmake")
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if (${STO_PREC_32} MATCHES "^[Tt]rue")
    add_compile_definitions(STO_PREC_32=1)
endif()
#add targets for libraries
add_executable( ${APP_EXE} src/cgs.cpp src/data_utils.cpp src/disp.cpp src/main.cpp )

#linking
target_link_libraries(${APP_EXE} PRIVATE Eigen3::Eigen meep png ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ctlgeom)
target_include_directories( ${APP_EXE} PRIVATE ${HDF5_INCLUDE_DIRS})
target_include_directories( ${APP_EXE} PRIVATE "~/install/include")

#check if we need to build tests
if(CMAKE_BUILD_TYPE MATCHES "^[Dd]ebug")
    #build and link the unit testing executable
    add_executable( ${TEST_EXE} src/cgs.cpp src/data_utils.cpp src/disp.cpp src/main_test.cpp )
    target_link_libraries(${TEST_EXE} PRIVATE Eigen3::Eigen meep ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
    target_include_directories( ${TEST_EXE} PRIVATE "~/install/include" )
    target_include_directories( ${TEST_EXE} PRIVATE "/usr/include/doctest" )
    target_include_directories( ${TEST_EXE} PRIVATE ${HDF5_INCLUDE_DIRS} )
    target_include_directories( ${TEST_EXE} PRIVATE "/usr/include/doctest")
    #copy files needed for unit testing

    add_custom_command(TARGET ${TEST_EXE} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
		   ${CMAKE_SOURCE_DIR}/tests $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests)
	       #configure_file(tests/test.geom tests/test.geom COPYONLY)
	       #configure_file(tests/test.conf tests/test.conf COPYONLY)
	       #configure_file(tests/run.geom tests/run.geom COPYONLY)
	       #configure_file(tests/run.conf tests/run.conf COPYONLY)
endif()

#allow the user to specify what type of geometry they'd like to use
if(NOT JUNCTION_TYPE)
    set(JUNCTION_TYPE Au_graphene_box)
endif()

message(CHECK_START "Using junction type: ${JUNCTION_TYPE}")
#copy the junction geometry files
configure_file(junctions/${JUNCTION_TYPE}/junc_template.geom junc_template.geom COPYONLY)
configure_file(junctions/${JUNCTION_TYPE}/params.conf params.conf COPYONLY)

#copy utility scripts that are used
configure_file(scripts/run.sh run.sh COPYONLY)
configure_file(scripts/utils.py utils.py COPYONLY)
configure_file(scripts/phases.py phases.py COPYONLY)
configure_file(scripts/set_dims.py set_dims.py COPYONLY)
configure_file(scripts/time_space.py time_space.py COPYONLY)
configure_file(scripts/check_enes.py check_enes.py COPYONLY)
