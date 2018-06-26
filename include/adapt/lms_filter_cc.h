/* -*- c++ -*- */
/*
 * Copyright 2018 <+YOU OR YOUR COMPANY+>.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#ifndef INCLUDED_ADAPT_LMS_FILTER_CC_H
#define INCLUDED_ADAPT_LMS_FILTER_CC_H

#include <adapt/api.h>
#include <gnuradio/sync_decimator.h>

namespace gr {
namespace adapt {

/*!
 * \brief Least Mean Squares Adaptive Filter (complex in/out)
 * \ingroup adapt
 *
 * \details
 * This block implements a complex LMS-based adaptive filter [1].
 *
 * [1] Widrow, Bernard, John McCool, and Michael Ball. "The complex LMS algorithm." Proceedings of
 * the IEEE 63.4 (1975): 719-720.
 */
class ADAPT_API lms_filter_cc : virtual public gr::sync_decimator {
    protected:
    virtual gr_complex error(const gr_complex& desired, const gr_complex& out) = 0;
    virtual void update_tap(gr_complex& tap, const gr_complex& in) = 0;

    public:
    typedef boost::shared_ptr<lms_filter_cc> sptr;

    /*!
     * Make an LMS adaptive filter
     *
     * \param first_input Specifies whether first input is reference or error signal (bool)
     * \param num_taps Number of taps in the filter (int)
     * \param mu Gain of the update loop (float)
     * \param skip Specifies how many samples are skipped between
     * successive filter updates (unsigned)
     * \param decimation Decimation rate of the filter (unsigned)
     * \param adapt Controls whether filter taps are being updated (bool)
     * \param reset Reset filter taps (bool)
     */
    static sptr make(bool first_input,
                     int num_taps,
                     float mu,
                     unsigned skip,
                     unsigned decimation,
                     bool adapt,
                     bool reset);

    virtual void set_taps(const std::vector<gr_complex>& taps) = 0;
    virtual const std::vector<gr_complex>& get_taps() = 0;
    virtual float get_mu() const = 0;
    virtual void set_mu(float mu) = 0;
    virtual unsigned get_skip() const = 0;
    virtual void set_skip(unsigned skip) = 0;
    virtual bool get_adapt() const = 0;
    virtual void set_adapt(bool adapt) = 0;
    virtual bool get_reset() const = 0;
    virtual void set_reset(bool reset) = 0;
};

} // namespace adapt
} // namespace gr

#endif /* INCLUDED_ADAPT_LMS_FILTER_CC_H */
