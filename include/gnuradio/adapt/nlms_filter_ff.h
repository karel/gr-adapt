/* -*- c++ -*- */
/*
 * Copyright 2018 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_ADAPT_NLMS_FILTER_FF_H
#define INCLUDED_ADAPT_NLMS_FILTER_FF_H

#include <gnuradio/adapt/api.h>
#include <gnuradio/sync_decimator.h>

namespace gr {
namespace adapt {

/*!
 * \brief Normalized Least Mean Squares Adaptive Filter (float in/out)
 * \ingroup adapt
 *
 * \details
 * This block implements an NLMS-based adaptive filter.
 */
class ADAPT_API nlms_filter_ff : virtual public gr::sync_decimator
{
protected:
    virtual float error(const float& desired, const float& out) = 0;
    virtual void update_tap(float& tap, const float& in) = 0;

public:
    typedef std::shared_ptr<nlms_filter_ff> sptr;

    /*!
     * Make a NLMS adaptive filter
     *
     * \param first_input Specifies whether first input is reference or error signal
     * (bool)
     * \param num_taps Number of taps in the filter (int)
     * \param mu Gain of the update loop (float)
     * \param skip Specifies how many samples are skipped between
     * successive filter updates (unsigned)
     * \param decimation Decimation rate of the filter (unsigned)
     * \param adapt Controls whether filter taps are being updated (bool)
     * \param bypass Bypass filter (bool)
     * \param reset Reset filter taps (bool)
     */
    static sptr make(bool first_input,
                     int num_taps,
                     float mu,
                     unsigned skip,
                     unsigned decimation,
                     bool adapt,
                     bool bypass,
                     bool reset);

    virtual void set_taps(const std::vector<float>& taps) = 0;
    virtual const std::vector<float>& get_taps() = 0;
    virtual float get_mu() const = 0;
    virtual void set_mu(float mu) = 0;
    virtual unsigned get_skip() const = 0;
    virtual void set_skip(unsigned skip) = 0;
    virtual bool get_adapt() const = 0;
    virtual void set_adapt(bool adapt) = 0;
    virtual bool get_bypass() const = 0;
    virtual void set_bypass(bool bypass) = 0;
    virtual bool get_reset() const = 0;
    virtual void set_reset(bool reset) = 0;
};

} // namespace adapt
} // namespace gr

#endif /* INCLUDED_ADAPT_NLMS_FILTER_FF_H */
