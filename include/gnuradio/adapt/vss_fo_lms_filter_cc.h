/* -*- c++ -*- */
/*
 * Copyright 2025 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_ADAPT_VSS_FO_LMS_FILTER_CC_H
#define INCLUDED_ADAPT_VSS_FO_LMS_FILTER_CC_H

#include <gnuradio/adapt/api.h>
#include <gnuradio/block.h>

namespace gr {
namespace adapt {

/*!
 * \brief Variable Step Sizes Frequency Offsets-Compensated Least Mean Squares
 * Adaptive filter (complex in/out)
 * \ingroup adapt
 *
 */
class ADAPT_API vss_fo_lms_filter_cc : virtual public gr::block
{
protected:
    virtual gr_complex error(const gr_complex& desired, const gr_complex& out) = 0;
    virtual void update_tap(gr_complex& tap, const gr_complex& in) = 0;

public:
    typedef std::shared_ptr<vss_fo_lms_filter_cc> sptr;

    /*!
     * Make a VSS-FO-LMS adaptive filter
     *
     * \param samp_rate Sampling rate (double)
     * \param num_taps Number of taps in the filter (int)
     * \param lambda_e Forgetting factor for estimation error averaging (float)
     * \param lambda_R Forgetting factor for cross-correlation vector averaging (float)
     * \param lambda_y Forgetting factor for input signal power averaging (float)
     * \param lambda_epsilon Forgetting factor for carrier frequency offset gradient
     * averaging (float)
     * \param lambda_eta Forgetting factor for sampling frequency offset gradient
     * averaging (float)
     * \param sigma2_v Variance of measurement noise (float)
     * \param step_size_limits Limits for variable step sizes [mu_cir_max, mu_cir_min,
     * mu_cfo_max, mu_cfo_min, mu_sfo_max, mu_sfo_min]
     * \param adapt Controls whether filter taps are being updated (bool)
     * \param reset Reset filter taps (bool).
     */
    static sptr make(double samp_rate,
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

    virtual float get_lambda_e() const = 0;
    virtual void set_lambda_e(float lambda_e) = 0;
    virtual float get_lambda_R() const = 0;
    virtual void set_lambda_R(float lambda_R) = 0;
    virtual float get_lambda_y() const = 0;
    virtual void set_lambda_y(float lambda_y) = 0;
    virtual float get_lambda_epsilon() const = 0;
    virtual void set_lambda_epsilon(float lambda_epsilon) = 0;
    virtual float get_lambda_eta() const = 0;
    virtual void set_lambda_eta(float lambda_eta) = 0;
    virtual float get_sigma2_v() const = 0;
    virtual void set_sigma2_v(float sigma2_v) = 0;
    virtual bool get_adapt() const = 0;
    virtual void set_adapt(bool adapt) = 0;
    virtual bool get_reset() const = 0;
    virtual void set_reset(bool reset) = 0;
};

} // namespace adapt
} // namespace gr

#endif /* INCLUDED_ADAPT_VSS_FO_LMS_FILTER_CC_H */
