find_package(PkgConfig)

PKG_CHECK_MODULES(PC_GR_ADAPT gnuradio-adapt)

FIND_PATH(
    GR_ADAPT_INCLUDE_DIRS
    NAMES gnuradio/adapt/api.h
    HINTS $ENV{ADAPT_DIR}/include
        ${PC_ADAPT_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    GR_ADAPT_LIBRARIES
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

include("${CMAKE_CURRENT_LIST_DIR}/gnuradio-adaptTarget.cmake")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GR_ADAPT DEFAULT_MSG GR_ADAPT_LIBRARIES GR_ADAPT_INCLUDE_DIRS)
MARK_AS_ADVANCED(GR_ADAPT_LIBRARIES GR_ADAPT_INCLUDE_DIRS)
