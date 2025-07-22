# ---------------------------------------------------------------------------
#  Tell CMake where Gurobi lives
# ---------------------------------------------------------------------------

# Add include + lib directories to the search path *early*
list(PREPEND CMAKE_PREFIX_PATH  "${GUROBI_HOME}")
list(PREPEND CMAKE_FIND_LIBRARY_PREFIXES "" "lib")
list(PREPEND CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

# ---------------------------------------------------------------------------
#  Look for the headers and the two libraries Gurobi needs
# ---------------------------------------------------------------------------
find_path(GUROBI_INCLUDE_DIR
        NAMES gurobi_c++.h
        HINTS "${GUROBI_HOME}/include")
find_library(GUROBI_CPP_LIB
        NAMES gurobi_c++
        HINTS "${GUROBI_HOME}/lib")
find_library(GUROBI_CORE_LIB
        NAMES gurobi gurobi120        # cover 9.0–11.*
        HINTS "${GUROBI_HOME}/lib")

if(NOT (GUROBI_INCLUDE_DIR AND GUROBI_CPP_LIB AND GUROBI_CORE_LIB))
    message(${GUROBI_INCLUDE_DIR}, ${GUROBI_CPP_LIB}, ${GUROBI_CORE_LIB})
    message(FATAL_ERROR "Could not find a complete Gurobi installation. "
            "Set GUROBI_HOME or pass -DGUROBI_HOME=… to cmake.")
endif()

# ---------------------------------------------------------------------------
#  Expose an IMPORTED target so the rest of the project can just do
#       target_link_libraries(foo PRIVATE Gurobi::gurobi)
# ---------------------------------------------------------------------------
add_library(Gurobi::gurobi INTERFACE IMPORTED)
set_target_properties(Gurobi::gurobi PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES      "${GUROBI_CPP_LIB};${GUROBI_CORE_LIB};pthread")
# ---------------------------------------------------------------------------
