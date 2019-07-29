INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_ADAPT adapt)

FIND_PATH(
    ADAPT_INCLUDE_DIRS
    NAMES adapt/api.h
    HINTS $ENV{ADAPT_DIR}/include
        ${PC_ADAPT_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    ADAPT_LIBRARIES
    NAMES gnuradio-adapt
    HINTS $ENV{ADAPT_DIR}/lib
        ${PC_ADAPT_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
          )

include("${CMAKE_CURRENT_LIST_DIR}/adaptTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ADAPT DEFAULT_MSG ADAPT_LIBRARIES ADAPT_INCLUDE_DIRS)
MARK_AS_ADVANCED(ADAPT_LIBRARIES ADAPT_INCLUDE_DIRS)
