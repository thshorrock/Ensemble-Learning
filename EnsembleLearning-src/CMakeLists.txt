## The Ensemble Learning library
project (EnsembleLearning)

##########################################################################
#  The Ensemble Learning library depends on the following boost packages:
#    - utility
#   A
##########################################################################

#set the libs to include
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

#set the include directory
INCLUDE_DIRECTORIES( "include" ${BOOST_INCLUDE_DIRS} )

#make the ensemble learning library
ADD_LIBRARY(EnsembleLearning  src/Random.cpp src/Builder.cpp src/Factory.cpp src/Placeholder.cpp )
target_link_libraries(EnsembleLearning  ${LIBS}) #link
install(DIRECTORY include/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -pg")
#set(CMAKE_BUILD_TYPE DEBUG)
install (TARGETS EnsembleLearning DESTINATION lib)
#Test
#add_subdirectory(test)


#Build Documentation
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${EnsembleLearningProject_SOURCE_DIR}/Doxyfile.in ${EnsembleLearningProject_SOURCE_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${EnsembleLearningProject_SOURCE_DIR}/Doxyfile
WORKING_DIRECTORY ${EnsembleLearningProject_SOURCE_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)