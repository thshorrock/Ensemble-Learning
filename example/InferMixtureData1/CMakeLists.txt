## Make the InferData example project
project (EnsembleLearning-InferMixture1)


#set the libs to include
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

#require the omp library
FIND_PACKAGE(OpenMP REQUIRED) 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

#include EL library
include_directories ("${EnsembleLearning_SOURCE_DIR}/include")
#include the includes of this project
include_directories ("${EnsembleLearning-InferMixture1_SOURCE_DIR}/../include")

set(CMAKE_BUILD_TYPE debug)

#build
add_executable(EnsembleLearning-InferMixture1  src/InferMixture1.cpp)
target_link_libraries (EnsembleLearning-InferMixture1 ${Boost_LIBRARIES})
target_link_libraries (EnsembleLearning-InferMixture1 ${boost_numeric_LIBRARY} ${boost_filesystem_LIBRARY} ${boost_program_options_LIBRARY} ${boost_system_LIBRARY}  EnsembleLearning)

install (TARGETS EnsembleLearning-InferMixture1 DESTINATION bin)