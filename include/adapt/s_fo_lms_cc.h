/* -*- c++ -*- */
/*
 * Copyright 2024 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_ADAPT_S_FO_LMS_CC_H
#define INCLUDED_ADAPT_S_FO_LMS_CC_H

#include <adapt/api.h>
#include <gnuradio/block.h>

namespace gr {
namespace adapt {

/*!
 * \brief Spline-Interpolated and Frequency Offsets-Compensated Least Mean Squares Adaptive filter
 * (complex in/out)
 * \ingroup adapt
 *
 */
class ADAPT_API s_fo_lms_cc : virtual public gr::block {
    protected:
    virtual gr_complex error(const gr_complex& desired, const gr_complex& out) = 0;
    virtual void update_tap(gr_complex& tap, const gr_complex& in) = 0;

    public:
    typedef std::shared_ptr<s_fo_lms_cc> sptr;

    /*!
     * Make an S-FO-LMS adaptive filter
     *
     * \param samp_rate Sampling rate (double)
     * \param num_taps Number of taps in the filter (int)
     * \param mu_cir Gain of the channel update loop (float)
     * \param mu_cfo Gain of the carrier frequency offset update loop (float)
     * \param mu_sfo Gain of the sampling frequency offset update loop (float)
     * \param mu_q Gain of the nonlinearities update loop (float)
     * \param adapt Controls whether filter taps are being updated (bool)
     * \param reset Reset filter taps (bool).
     */
    static sptr make(double samp_rate,
                     int num_taps,
                     float mu_cir,
                     float mu_cfo,
                     float mu_sfo,
                     float mu_q,
                     bool adapt,
                     bool reset);

    virtual float get_mu_cir() const = 0;
    virtual void set_mu_cir(float mu_cir) = 0;
    virtual float get_mu_cfo() const = 0;
    virtual void set_mu_cfo(float mu_cfo) = 0;
    virtual float get_mu_sfo() const = 0;
    virtual void set_mu_sfo(float mu_sfo) = 0;
    virtual float get_mu_q() const = 0;
    virtual void set_mu_q(float mu_q) = 0;
    virtual bool get_adapt() const = 0;
    virtual void set_adapt(bool adapt) = 0;
    virtual bool get_reset() const = 0;
    virtual void set_reset(bool reset) = 0;
};

} // namespace adapt
} // namespace gr

#endif /* INCLUDED_ADAPT_S_FO_LMS_CC_H */
