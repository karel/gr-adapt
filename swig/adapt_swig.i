/* -*- c++ -*- */

#define ADAPT_API

%include "gnuradio.i"           // the common stuff

//load generated python docstrings
%include "adapt_swig_doc.i"

%{
#include "adapt/lms_filter_ff.h"
#include "adapt/lms_filter_cc.h"
#include "adapt/nlms_filter_ff.h"
#include "adapt/nlms_filter_cc.h"
#include "adapt/rls_filter_ff.h"
#include "adapt/rls_filter_cc.h"
#include "adapt/qrd_rls_filter_ff.h"
#include "adapt/qrd_rls_filter_cc.h"
#include "adapt/iqrd_rls_filter_ff.h"
#include "adapt/iqrd_rls_filter_cc.h"
%}

%include "adapt/lms_filter_ff.h"
GR_SWIG_BLOCK_MAGIC2(adapt, lms_filter_ff);
%include "adapt/lms_filter_cc.h"
GR_SWIG_BLOCK_MAGIC2(adapt, lms_filter_cc);
%include "adapt/nlms_filter_ff.h"
GR_SWIG_BLOCK_MAGIC2(adapt, nlms_filter_ff);
%include "adapt/nlms_filter_cc.h"
GR_SWIG_BLOCK_MAGIC2(adapt, nlms_filter_cc);
%include "adapt/rls_filter_ff.h"
GR_SWIG_BLOCK_MAGIC2(adapt, rls_filter_ff);
%include "adapt/rls_filter_cc.h"
GR_SWIG_BLOCK_MAGIC2(adapt, rls_filter_cc);
%include "adapt/qrd_rls_filter_ff.h"
GR_SWIG_BLOCK_MAGIC2(adapt, qrd_rls_filter_ff);
%include "adapt/qrd_rls_filter_cc.h"
GR_SWIG_BLOCK_MAGIC2(adapt, qrd_rls_filter_cc);
%include "adapt/iqrd_rls_filter_ff.h"
GR_SWIG_BLOCK_MAGIC2(adapt, iqrd_rls_filter_ff);
%include "adapt/iqrd_rls_filter_cc.h"
GR_SWIG_BLOCK_MAGIC2(adapt, iqrd_rls_filter_cc);
