## Make the ICA example project
project (ICA)


#set the libs to include
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

#require the omp library
FIND_PACKAGE(OpenMP REQUIRED) 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -pg")
#set(CMAKE_BUILD_TYPE DEBUG)

#include EL library
include_directories ("${EnsembleLearning_SOURCE_DIR}/include")
#include the includes of this project
include_directories ("${ICA_SOURCE_DIR}/../include")

#build
add_executable(ICA  ../src/ICA.cpp)
target_link_libraries (ICA ${Boost_LIBRARIES})
#set_target_properties (ICA PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
#set_target_properties (ICA PROPERTIES LINKER_FLAGS ${OpenMP_C_FLAGS})
target_link_libraries (ICA ${boost_numeric_LIBRARY} ${boost_filesystem_LIBRARY} ${boost_program_options_LIBRARY} ${boost_system_LIBRARY}  EnsembleLearning)
#boost_report_pretty("project src" ${EnsembleLearning_SOURCE_DIR} )

install (TARGETS ICA DESTINATION bin)