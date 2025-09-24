/* -*- c++ -*- */
/*
 * Copyright 2023 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_ADAPT_FO_LMS_FILTER_CC_IMPL_H
#define INCLUDED_ADAPT_FO_LMS_FILTER_CC_IMPL_H

#include <gnuradio/adapt/fo_lms_filter_cc.h>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/filter/mmse_fir_interpolator_cc.h>
#include <gnuradio/gr_complex.h>

namespace gr {
namespace adapt {

class fo_lms_filter_cc_impl : public fo_lms_filter_cc, filter::kernel::fir_filter_ccc
{
private:
    std::vector<gr_complex> d_new_taps;
    bool d_updated = false;
    const gr::filter::mmse_fir_interpolator_cc d_interpolator;
    double d_samp_rate;
    gr_complex d_error, d_out_previous;
    gr_complex* d_y;
    float d_mu_cir, d_mu_cfo, d_mu_sfo;
    float d_cfo = 0.0, d_sfo = 0.0;
    float d_p = 0.0, d_t = 0.0;
    int d_M;
    bool d_adapt, d_reset;

protected:
    gr_complex error(const gr_complex& desired, const gr_complex& out) override;
    void update_tap(gr_complex& tap, const gr_complex& in) override;

public:
    fo_lms_filter_cc_impl(double d_samp_rate,
                          int num_taps,
                          float mu_cir,
                          float mu_cfo,
                          float mu_sfo,
                          bool adapt,
                          bool reset);
    ~fo_lms_filter_cc_impl();

    float get_mu_cir() const override;
    void set_mu_cir(float mu_cir) override;
    float get_mu_cfo() const override;
    void set_mu_cfo(float mu_cfo) override;
    float get_mu_sfo() const override;
    void set_mu_sfo(float mu_sfo) override;
    bool get_adapt() const override;
    void set_adapt(bool adapt) override;
    bool get_reset() const override;
    void set_reset(bool reset) override;

    // Where all the action really happens
    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items) override;
};

} // namespace adapt
} // namespace gr

#endif /* INCLUDED_ADAPT_FO_LMS_FILTER_CC_IMPL_H */
