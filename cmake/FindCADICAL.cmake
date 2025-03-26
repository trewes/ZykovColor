# This module finds the cadical header and libraries.
#
# User can give CADICAL_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  CADICAL_FOUND              - Set to false, or undefined, if cadical isn't found.
#  CADICAL_INCLUDE_DIR        - include directory
#  CADICAL_LIBRARY            - library files

# use given hint directory or look in parent/cadical folder
IF(NOT DEFINED CADICAL_ROOT_DIR)
    SET(CADICAL_ROOT_DIR ${CMAKE_SOURCE_DIR}/../cadical)
ENDIF()

FIND_PATH(CADICAL_INCLUDE_DIR cadical.hpp
        PATH_SUFFIXES src
        PATHS ${CADICAL_ROOT_DIR}
)
MESSAGE(STATUS "CADICAL Include Dir: ${CADICAL_INCLUDE_DIR}")


FIND_LIBRARY(CADICAL_LIBRARY
        NAMES cadical libcadical.a
        PATH_SUFFIXES lib build
        PATHS ${CADICAL_ROOT_DIR}
)
MESSAGE(STATUS "CADICAL Library: ${CADICAL_LIBRARY}")


# check whether required things have been found and set CADICAL_FOUND accordingly
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CADICAL DEFAULT_MSG
        CADICAL_INCLUDE_DIR CADICAL_LIBRARY)

