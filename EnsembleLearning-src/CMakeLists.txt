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
install (TARGETS EnsembleLearning DESTINATION lib) #install
#install the files
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
install (FILES ${files}  DESTINATION include)
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/node/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/node )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/node/factor/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/node/factor )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/node/variable/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/node/variable )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/message/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/message )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/exponential_model/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/exponential_model )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/calculation_tree/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/calculation_tree )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/detail/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/detail )
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/include/EnsembleLearning/exception/*.hpp")
install (FILES ${files}  DESTINATION include/EnsembleLearning/exception/)

#Test
add_subdirectory(test)