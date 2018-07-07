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

#ifndef INCLUDED_ADAPT_NLMS_FILTER_CC_IMPL_H
#define INCLUDED_ADAPT_NLMS_FILTER_CC_IMPL_H

#include "config.h"
#ifdef ARMADILLO_FOUND
#include <armadillo>
#endif
#include <adapt/nlms_filter_cc.h>
#include <gnuradio/filter/fir_filter.h>

namespace gr {
namespace adapt {

class nlms_filter_cc_impl : public nlms_filter_cc, filter::kernel::fir_filter_ccc {
    private:
#ifdef ARMADILLO_FOUND
    arma::cx_fvec d_taps;
    arma::cx_fvec d_new_taps;
#else
    std::vector<gr_complex> d_new_taps;
#endif // ARMADILLO_FOUND
    bool d_first_input;
    bool d_updated;
    gr_complex d_error, d_power;
    float d_mu, d_epsilon;
    unsigned d_skip, d_i;
    bool d_adapt, d_reset;

    protected:
    gr_complex error(const gr_complex& desired, const gr_complex& out) override;
    void update_tap(gr_complex& tap, const gr_complex& in) override;

    public:
    nlms_filter_cc_impl(bool first_input,
                        int num_taps,
                        float mu,
                        unsigned skip,
                        unsigned decimation,
                        bool adapt,
                        bool reset);
    ~nlms_filter_cc_impl() = default;

    const std::vector<gr_complex>& get_taps() override;
    void set_taps(const std::vector<gr_complex>& new_taps) override;
    float get_mu() const override;
    void set_mu(float mu) override;
    unsigned get_skip() const override;
    void set_skip(unsigned skip) override;
    bool get_adapt() const override;
    void set_adapt(bool adapt) override;
    bool get_reset() const override;
    void set_reset(bool reset) override;

    // Where all the action really happens
    int work(int noutput_items,
             gr_vector_const_void_star& input_items,
             gr_vector_void_star& output_items) override;
};

} // namespace adapt
} // namespace gr

#endif /* INCLUDED_ADAPT_NLMS_FILTER_CC_IMPL_H */