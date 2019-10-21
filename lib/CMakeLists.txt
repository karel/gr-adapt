# Copyright 2011,2012,2016 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

########################################################################
# Setup library
########################################################################
include(GrPlatform) #define LIB_SUFFIX

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
list(APPEND adapt_sources
	lms_filter_ff_impl.cc
    lms_filter_cc_impl.cc
	nlms_filter_ff_impl.cc
    nlms_filter_cc_impl.cc
    rls_filter_ff_impl.cc
    rls_filter_cc_impl.cc
    qrd_rls_filter_ff_impl.cc
    qrd_rls_filter_cc_impl.cc
    iqrd_rls_filter_ff_impl.cc
    iqrd_rls_filter_cc_impl.cc
    )

set(adapt_sources "${adapt_sources}" PARENT_SCOPE)
if(NOT adapt_sources)
	MESSAGE(STATUS "No C++ sources... skipping lib/")
	return()
endif(NOT adapt_sources)

add_library(gnuradio-adapt SHARED ${adapt_sources})
target_link_libraries(gnuradio-adapt ${VOLK_LIBRARIES} ${Boost_LIBRARIES} ${GNURADIO_ALL_LIBRARIES})
if(LAPACK_FOUND AND ARMADILLO_FOUND)
	target_link_libraries(gnuradio-adapt ${LAPACK_LIBRARIES} ${ARMADILLO_LIBRARIES})
endif(LAPACK_FOUND AND ARMADILLO_FOUND)
set_target_properties(gnuradio-adapt PROPERTIES DEFINE_SYMBOL "gnuradio_adapt_EXPORTS")

if(APPLE)
    set_target_properties(gnuradio-adapt PROPERTIES
        INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib"
    )
endif(APPLE)

########################################################################
# Install built library files
########################################################################
include(GrMiscUtils)
GR_LIBRARY_FOO(gnuradio-adapt RUNTIME_COMPONENT "adapt_runtime" DEVEL_COMPONENT "adapt_devel")

########################################################################
# Build and register unit test
########################################################################
include(GrTest)

include_directories(${CPPUNIT_INCLUDE_DIRS})

list(APPEND test_adapt_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_adapt.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qa_adapt.cc
)

add_executable(test-adapt ${test_adapt_sources})

target_link_libraries(
  test-adapt
  ${GNURADIO_RUNTIME_LIBRARIES}
  ${Boost_LIBRARIES}
  ${CPPUNIT_LIBRARIES}
  gnuradio-adapt
)

GR_ADD_TEST(test_adapt test-adapt)

########################################################################
# Print summary
########################################################################
message(STATUS "Using install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Building for version: ${VERSION} / ${LIBVER}")

