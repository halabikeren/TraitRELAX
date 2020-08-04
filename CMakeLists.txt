# CMake script for TraitRELAX
# Author:
#   Keren Halabi
# Created: 02/08/2020

cmake_minimum_required (VERSION 2.8.11)
project (traitrelax CXX)

# Compile options
set (CMAKE_CXX_FLAGS "-std=c++11 -Wall -Weffc++ -Wshadow -Wconversion")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NO_DEP_CHECK)
  MESSAGE("-- Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(${PROJECT_NAME}_VERSION_CURRENT "12")
SET(${PROJECT_NAME}_VERSION_REVISION "0")
SET(${PROJECT_NAME}_VERSION_AGE "0")

# Effective version number computation
MATH(EXPR ${PROJECT_NAME}_VERSION_MAJOR "${${PROJECT_NAME}_VERSION_CURRENT} - ${${PROJECT_NAME}_VERSION_AGE}")
SET(${PROJECT_NAME}_VERSION_MINOR ${${PROJECT_NAME}_VERSION_AGE})
SET(${PROJECT_NAME}_VERSION_PATCH ${${PROJECT_NAME}_VERSION_REVISION})
SET(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

set (PROJECT_VERSION ${${PROJECT_NAME}_VERSION})

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-phyl 12.0.0 REQUIRED)

# CMake package
set (cmake-package-location ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
include (CMakePackageConfigHelpers)
configure_package_config_file (
  package.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/package.cmake
  INSTALL_DESTINATION ${cmake-package-location}
  )
write_basic_package_version_file (
  ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config.cmake)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config-version.cmake)

# Define the libraries
add_subdirectory (lib)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "traitrelax")
SET(CPACK_PACKAGE_VERSION "1.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "TraitRELAX")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES 
       "/build/" 
       "/html/" 
       "/TraitRELAX.tag" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "lib${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

IF(NOT NO_DEP_CHECK)
IF (UNIX)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist)

ENDIF()

SET(CTEST_UPDATE_TYPE git)
SET(UPDATE_COMMAND "git")
SET(UPDATE_OPTIONS "")

enable_testing ()
include (CTest)
if (BUILD_TESTING)
  add_subdirectory (TraitRELAX)
endif (BUILD_TESTING)

ENDIF(NOT NO_DEP_CHECK)