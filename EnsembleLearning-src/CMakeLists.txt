## The Ensemble Learning library
project (EnsembleLearning)

##########################################################################
#  The Ensemble Learning library depends on the following boost packages:
#    - utility
#    - assign
#   A
##########################################################################

#set the libs to include
set(LIBS ${LIBS} ${GSL_LIBRARIES} ${GSLCBLAS_LIBRARIES})

#set the include directory
INCLUDE_DIRECTORIES( "include" ${BOOST_INCLUDE_DIRS} )

#make the ensemble learning library
ADD_LIBRARY(EnsembleLearning  src/Builder.cpp)
target_link_libraries(EnsembleLearning  ${LIBS}) #link
install(DIRECTORY include/ DESTINATION include
          FILES_MATCHING PATTERN "*.hpp")

#Test
add_subdirectory(test)


#Build Documentation
# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)