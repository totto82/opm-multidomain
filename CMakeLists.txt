# -*- mode: cmake; tab-width: 2; indent-tabs-mode: nil; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

# Mandatory call to project
project(opm-multidomain C CXX)

cmake_minimum_required (VERSION 2.8)

# add the current projects cmake module directory to the search
# path. This is not required anymore once support for federated builds
# is merged.
list(INSERT CMAKE_MODULE_PATH 0 "${PROJECT_SOURCE_DIR}/cmake/Modules")
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)

# Testing
include(OpmSatellites)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include(${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include(CMakeLists_files.cmake)

macro (config_hook)
  opm_need_version_of("dune-common")
  opm_need_version_of("dune-geometry")
  opm_need_version_of("dune-grid")
  opm_need_version_of("dune-localfunctions")
  opm_need_version_of("dune-alugrid")
  opm_need_version_of("dune-istl")
  opm_need_version_of("dune-fem")
endmacro (config_hook)

macro (files_hook)
endmacro (files_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (tests_hook)
endmacro (tests_hook)

# For eWoms, we set the BUILD_TESTING variable to off by default. The
# reason is that builing the tests takes quite some time and they are
# not very interesting for people who just want to use eWoms but not
# develop it. To enable testing, either pass --enable-testing to the
# configure script or call cmake with -DBUILD_TESTING=ON
set(MY_BUILD_TESTING "${BUILD_TESTING}")
set(BUILD_TESTING "ON" CACHE BOOL "Build the tests" FORCE)

# all setup common to the OPM library modules
include (OpmLibMain)

# restore the value of the BUILD_TESTING variable to its original
# value. This is required because OpmLibMain thinks it is smarter
# appropriate and sets it...
if (NOT MY_BUILD_TESTING)
  set(BUILD_TESTING "${MY_BUILD_TESTING}" CACHE BOOL "Build the tests" FORCE)
endif()

# copy the data files (i.e., grid files, etc) needed to run the unit
# tests from the source to the build directory
file(GLOB datafiles "tests/data/*")
foreach(file ${datafiles})
  file(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
endforeach()


opm_set_test_driver("${PROJECT_SOURCE_DIR}/bin/runtest.sh" "--simulation")
opm_set_test_default_working_directory("${PROJECT_BINARY_DIR}")

set(CMAKE_PROJECT_NAME "${PROJECT_NAME}")
add_custom_target(opm-models_prepare)

# add targets for all tests of the models.
foreach(tapp 
    )
  opm_add_test(${tapp})
endforeach()
