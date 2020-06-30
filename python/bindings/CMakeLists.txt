# Copyright 2020 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

GR_PYTHON_CHECK_MODULE_RAW(
    "pygccxml"
    "import pygccxml"
    PYGCCXML_FOUND
    )

include(GrPybind)

########################################################################
# Python Bindings
########################################################################
list(APPEND adapt_python_files
    lms_filter_ff_python.cc
    lms_filter_cc_python.cc
    nlms_filter_ff_python.cc
    nlms_filter_cc_python.cc
    qrd_rls_filter_ff_python.cc
    qrd_rls_filter_cc_python.cc
    rls_filter_ff_python.cc
    rls_filter_cc_python.cc
    iqrd_rls_filter_ff_python.cc
    iqrd_rls_filter_cc_python.cc python_bindings.cc)

GR_PYBIND_MAKE_OOT(adapt 
   ../..
   gr::adapt
   "${adapt_python_files}")

install(TARGETS adapt_python DESTINATION ${GR_PYTHON_DIR}/adapt COMPONENT pythonapi)
