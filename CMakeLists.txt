
cmake_minimum_required (VERSION 2.6)
project (EnsembleLearningProject)

# The version number.
set (EnsembleLearning_VERSION_MAJOR 0)
set (EnsembleLearning_VERSION_MINOR 1)
set (EnsembleLearning_PATCH_LEVEL 0 )

# project options
OPTION( BUILD_SHARED_LIBS "Set to OFF to build static libraries" ON )
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )
OPTION( BUILD_EXAMPLES "Set to OFF to skip building the examples" ON )
OPTION( BUILD_MISSING_DEPENDANCIES "Set to OFF to skip building missing dependencies (you might want to build the latest versions yourself)" ON )

# Put the libaries and binaries that get built into directories at the
# top of the build tree rather than in hard-to-find leaf
# directories.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


##########################################################################
# The libraries and executables depend upon boost and gsl library.
# If these already exist on system then want to use the existing install.
#  If not (and if BUILD_MISSING_DEPENDANCIES is ON) then want to build them.
###########################################################################

##########################################################################
#  BOOST FIRST;
##########################################################################
FIND_PACKAGE(Boost COMPONENTS system program_options test_exec_monitor filesystem REQUIRED) # do we have any boost at all?

IF(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
ELSE()
  IF(BUILD_MISSING_DEPENDANCIES)
  #At present the following is redundent (GSL required)
  #  This is because it is not working yet, when it is done the required will be removed
    message(STATUS "BUILDING call_traits library")  
    add_subdirectory(boost)
    include_directories(boost) 
    set(Boost_filesystem_LIBRARY ${CMAKE_BINARY_DIR}/lib/boost_filesytem-mt)
    set(Boost_program_options_LIBRARY ${CMAKE_BINARY_DIR}/lib/boost_program_options-mt)
  ELSE()
    message(STATUS "Could not find boost.  If you have boost you need to tell me where it is by using the ... command.\n  Alternatively, you may: \n\t - build the latest boost libraries yourself (Recommended),\n\t - Or let me build a (possibly old) version myself.\n\n  For this final option, you need the BUILD_MISSING_DEPENDANCIES option to be ON.")
  ENDIF()
ENDIF()

##########################################################################
#  Now gsl subset
##########################################################################

find_package(GSL REQUIRED)
IF("${GSL_FOUND}")
  include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})
ELSE("${GSL_FOUND}")  
  #At present the following is redundent (GSL required)
  #  This is because it is not working yet, when it is done the required will be removed
  MESSAGE("GSL NOT FOUND -- BUILDING FROM SOURCE")
  add_subdirectory(gsl)
ENDIF("${GSL_FOUND}" )

##########################################################################
# Ensemble Learning modules               
# Do in the same style as Boost.CMake so that uniform style if libraries not found      #
##########################################################################

include(boost/tools/build/CMake/BoostUtils.cmake)
message(STATUS "")
colormsg(_HIBLUE_ "Ensemble Learning Project, starting ...")
message(STATUS "")

##########################################################################
#  Build the EL library 
##########################################################################

enable_testing()

add_subdirectory(EnsembleLearning-src)

IF( BUILD_EXAMPLES )
  

  message(STATUS "")
  colormsg(_HIBLUE_ "Example starting")
  message(STATUS "")



  add_subdirectory(example/InferData)
  add_subdirectory(example/ICA/build)
  

ENDIF( BUILD_EXAMPLES )


INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Ensemble Learning Project")
SET(CPACK_PACKAGE_VENDOR "T.H. Shorrock")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE_1_0.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  #SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\ICA.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} EnsembleLearning")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\git@github.com:thshorrock/Ensemble-Learning.git")
  #SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.my-personal-home-page.com")
  SET(CPACK_NSIS_CONTACT "t.h.shorrock@gmail.com")
  SET(CPACK_NSIS_MODIFY_PATH ON)
ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_STRIP_FILES "bin/ICA")
  SET(CPACK_SOURCE_STRIP_FILES "")
ENDIF(WIN32 AND NOT UNIX)
SET(CPACK_PACKAGE_EXECUTABLES "ICA" "EnsembleLearning")
INCLUDE(CPack)