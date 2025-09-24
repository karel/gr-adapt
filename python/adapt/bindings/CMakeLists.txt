# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

########################################################################
# Check if there is C++ code at all
########################################################################
if(NOT adapt_sources)
    message(STATUS "No C++ sources... skipping python bindings")
    return()
endif(NOT adapt_sources)

########################################################################
# Check for pygccxml
########################################################################
gr_python_check_module_raw("pygccxml" "import pygccxml" PYGCCXML_FOUND)

include(GrPybind)

########################################################################
# Python Bindings
########################################################################

list(APPEND adapt_python_files
    lms_filter_ff_python.cc
    lms_filter_cc_python.cc
    nlms_filter_ff_python.cc
    nlms_filter_cc_python.cc
    rls_filter_ff_python.cc
    rls_filter_cc_python.cc
    qrd_rls_filter_ff_python.cc
    qrd_rls_filter_cc_python.cc
    iqrd_rls_filter_ff_python.cc
    iqrd_rls_filter_cc_python.cc
    fo_lms_filter_cc_python.cc
    s_fo_lms_filter_cc_python.cc python_bindings.cc)

gr_pybind_make_oot(adapt ../../.. gr::adapt "${adapt_python_files}")

# copy bindings extension for use in QA test module
add_custom_command(
    TARGET adapt_python
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:adapt_python>
            ${PROJECT_BINARY_DIR}/test_modules/gnuradio/adapt/)

install(
    TARGETS adapt_python
    DESTINATION ${GR_PYTHON_DIR}/gnuradio/adapt
    COMPONENT pythonapi)
