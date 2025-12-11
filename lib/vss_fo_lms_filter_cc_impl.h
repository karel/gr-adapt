/* -*- c++ -*- */
/*
 * Copyright 2025 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_ADAPT_VSS_FO_LMS_FILTER_CC_IMPL_H
#define INCLUDED_ADAPT_VSS_FO_LMS_FILTER_CC_IMPL_H

#include <gnuradio/adapt/vss_fo_lms_filter_cc.h>
#include <gnuradio/filter/fir_filter.h>
#include <gnuradio/filter/mmse_fir_interpolator_cc.h>
#include <gnuradio/gr_complex.h>
#include <vector>

namespace gr {
namespace adapt {

class vss_fo_lms_filter_cc_impl : public vss_fo_lms_filter_cc,
                                  filter::kernel::fir_filter_ccc
{
private:
    std::vector<gr_complex> d_new_taps;
    bool d_updated = false;
    const gr::filter::mmse_fir_interpolator_cc d_interpolator;
    double d_samp_rate;
    gr_complex d_error, d_out_previous;
    gr_complex* d_y;
    float d_mu_cir = 0.0, d_mu_cfo = 0.0, d_mu_sfo = 0.0;
    float d_mu_cir_min, d_mu_cir_max, d_mu_cfo_min, d_mu_cfo_max, d_mu_sfo_min,
        d_mu_sfo_max;
    std::vector<float> d_mu_cfos, d_mu_sfos;
    int d_M;
    float d_lambda_e, d_lambda_R, d_lambda_y, d_lambda_epsilon, d_lambda_eta;
    float d_sigma2_e = 1.0, d_sigma2_y = 1.0, d_sigma2_v = 0.0;
    gr_complex* d_Rye;
    float d_g_cfo = 0.0, d_g_sfo = 0.0;
    float d_cfo = 0.0, d_sfo = 0.0;
    float d_p = 0.0, d_t = 0.0;
    bool d_adapt, d_reset;

protected:
    gr_complex error(const gr_complex& desired, const gr_complex& out) override;
    void update_tap(gr_complex& tap, const gr_complex& in) override;

public:
    vss_fo_lms_filter_cc_impl(double samp_rate,
                              int num_taps,
                              float lambda_e,
                              float lambda_R,
                              float lambda_y,
                              float lambda_epsilon,
                              float lambda_eta,
                              float sigma2_v,
                              const std::vector<float>& step_size_limits,
                              bool adapt,
                              bool reset);
    ~vss_fo_lms_filter_cc_impl();

    float get_lambda_e() const override;
    void set_lambda_e(float lambda_e) override;
    float get_lambda_R() const override;
    void set_lambda_R(float lambda_R) override;
    float get_lambda_y() const override;
    void set_lambda_y(float lambda_y) override;
    float get_lambda_epsilon() const override;
    void set_lambda_epsilon(float lambda_epsilon) override;
    float get_lambda_eta() const override;
    void set_lambda_eta(float lambda_eta) override;
    float get_sigma2_v() const override;
    void set_sigma2_v(float sigma2_v) override;
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

#endif /* INCLUDED_ADAPT_VSS_FO_LMS_FILTER_CC_IMPL_H */
