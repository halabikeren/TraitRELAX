# CMake script for TraitRELAX program, based on the script from TestNH (https://github.com/BioPP/testnh)

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
  MESSAGE(STATUS "Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

# static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE(STATUS "Static linkage requested.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
ENDIF()


# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-phyl 12.0.0 REQUIRED)

# Subdirectories
add_subdirectory (TraitRELAX)

ENDIF(NO_DEP_CHECK)
