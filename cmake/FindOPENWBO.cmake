# This module finds the Open-WBO header and libraries.
#
# User can give OPENWBO_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  OPENWBO_FOUND              - Set to false, or undefined, if open-wbo isn't found.
#  OPENWBO_INCLUDE_DIRS       - include directories
#  OPENWBO_LIBRARY            - library files

# use given hint directory or look in parent/open-wbo folder
IF(NOT DEFINED OPENWBO_ROOT_DIR)
    SET(OPENWBO_ROOT_DIR ${CMAKE_SOURCE_DIR}/../open-wbo)
ENDIF()

# For this project we need the encoder header and the glucose headers
FIND_PATH(OPENWBO_ENCODER_DIR Encoder.h
        PATHS ${OPENWBO_ROOT_DIR})
FIND_PATH(OPENWBO_SOLVERS_DIR
        NAMES glucose4.2
        PATH_SUFFIXES solvers
        HINTS ${OPENWBO_ROOT_DIR}/solvers)
SET(OPENWBO_GLUCOSE_DIR ${OPENWBO_SOLVERS_DIR}/glucose4.2/)

# Combine the found directories into a single variable
SET(OPENWBO_INCLUDE_DIRS  ${OPENWBO_ENCODER_DIR} ${OPENWBO_GLUCOSE_DIR})

# open-wbo compiles with the macro NSPACE depending on which solver is used, have to repeat it here too
ADD_DEFINITIONS(-DNSPACE=Glucose) #TODO nicer way to do this?

MESSAGE(STATUS "OPENWBO Include Dirs: ${OPENWBO_INCLUDE_DIRS}")

FIND_LIBRARY(OPENWBO_LIBRARY
        NAMES lib lib.a
        PATH_SUFFIXES
        PATHS ${OPENWBO_ROOT_DIR}
)
MESSAGE(STATUS "OPENWBO Library: ${OPENWBO_LIBRARY}")


# check whether required things have been found and set OPENWBO_FOUND accordingly
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OPENWBO DEFAULT_MSG
        OPENWBO_INCLUDE_DIRS OPENWBO_ENCODER_DIR OPENWBO_GLUCOSE_DIR OPENWBO_LIBRARY)

