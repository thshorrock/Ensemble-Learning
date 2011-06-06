## Make the InferData example project
project (EnsembleLearning-InferData1)


#set the libs to include
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

#require the omp library
FIND_PACKAGE(OpenMP REQUIRED) 
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

#include EL library
include_directories ("${EnsembleLearning_SOURCE_DIR}/include")
#include the includes of this project
include_directories ("${EnsembleLearning-InferData1_SOURCE_DIR}/../include")

set(CMAKE_BUILD_TYPE debug)

#build
add_executable(EnsembleLearning-InferData1  src/InferData1.cpp)
target_link_libraries (EnsembleLearning-InferData1 ${Boost_LIBRARIES})
target_link_libraries (EnsembleLearning-InferData1 ${boost_numeric_LIBRARY} ${boost_filesystem_LIBRARY} ${boost_program_options_LIBRARY} ${boost_system_LIBRARY}  EnsembleLearning)

install (TARGETS EnsembleLearning-InferData1 DESTINATION bin)