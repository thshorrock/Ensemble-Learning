

##############################
#  Make the test binary
#############################


project(test_EnsembleLearning)


FIND_PACKAGE(OpenMP REQUIRED) 
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_BUILD_TYPE DEBUG)

IF(BUILD_COVERAGE)
 IF (CMAKE_COMPILER_IS_GNUCC)
   set(CMAKE_BUILD_TYPE DEBUG)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
   set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
  
    ADD_CUSTOM_COMMAND(TARGET helloworld PRE_BUILD
        COMMAND ${CMAKE_SOURCE_DIR}/../scripts/test4ccache ${CMAKE_C_COMPILER}
        COMMENT "Testing for ccache"
        )
    SET(CMAKE_C_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage --coverage ")
    SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fprofile-arcs -ftest-coverage --coverage ")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage  --coverage -lgcov")
    # add lcov target
    include(${CMAKE_SOURCE_DIR}/cmake/lcov.cmake)
    #ADD_DEPENDENCIES(lcov helloworld)
 ENDIF()
#require the omp library
ENDIF(BUILD_COVERAGE)

add_executable(test_EnsembleLearning test_EnsembleLearning.cpp)
target_link_libraries (test_EnsembleLearning EnsembleLearning)
target_link_libraries (test_EnsembleLearning ${Boost_LIBRARIES})
target_link_libraries (test_EnsembleLearning ${boost_numeric_LIBRARY} ${boost_filesystem_LIBRARY} ${boost_program_options_LIBRARY} ${boost_system_LIBRARY} ${boost_test_LIBRARY}  EnsembleLearning  ${boost_test_exec_monitor_LIBRARY} )
add_test(test-EnsembleLearning-lib ${CMAKE_BINARY_DIR}/bin/test_EnsembleLearning Hello )