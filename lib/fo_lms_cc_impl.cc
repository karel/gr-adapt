/* -*- c++ -*- */
/*
 * Copyright 2023 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */


#include "fo_lms_cc_impl.h"
#include <cmath>
#include <gnuradio/block.h>
#include <gnuradio/gr_complex.h>
#include <gnuradio/io_signature.h>
#include <gnuradio/logger.h>
#include <string>
#include <volk/volk.h>

namespace gr {
namespace adapt {

using namespace filter::kernel;

using input_type = gr_complex;
using output_type = gr_complex;
fo_lms_cc::sptr
fo_lms_cc::make(double samp_rate, int num_taps, float mu_cir, float mu_cfo, float mu_sfo) {
    return gnuradio::make_block_sptr<fo_lms_cc_impl>(samp_rate, num_taps, mu_cir, mu_cfo, mu_sfo);
}


/*
 * The private constructor
 */
fo_lms_cc_impl::fo_lms_cc_impl(
    double samp_rate, int num_taps, float mu_cir, float mu_cfo, float mu_sfo)
    : gr::block("fo_lms_cc",
                gr::io_signature::make(2, 2, sizeof(input_type)),
                gr::io_signature::makev(
                    1,
                    4,
                    std::vector<int>{
                        sizeof(output_type), sizeof(output_type), sizeof(float), sizeof(float)})),
      fir_filter_ccc(std::vector<gr_complex>(num_taps, gr_complex(0, 0))), d_samp_rate(samp_rate),
      d_mu_cir(mu_cir), d_mu_cfo(mu_cfo), d_mu_sfo(mu_sfo) {
    const int alignment_multiple = volk_get_alignment() / sizeof(input_type);
    set_alignment(std::max(1, alignment_multiple));

    d_y_0 = new gr_complex[num_taps]{};
    d_y_1 = new gr_complex[num_taps]{};
    d_y_2 = new gr_complex[num_taps]{};

    set_history(num_taps);
}

fo_lms_cc_impl::~fo_lms_cc_impl() {
    delete[] d_y_0;
    delete[] d_y_1;
    delete[] d_y_2;
}

float fo_lms_cc_impl::get_mu_cir() const { return d_mu_cir; }

void fo_lms_cc_impl::set_mu_cir(float mu_cir) { d_mu_cir = mu_cir; }

float fo_lms_cc_impl::get_mu_cfo() const { return d_mu_cfo; }

void fo_lms_cc_impl::set_mu_cfo(float mu_cfo) { d_mu_cfo = mu_cfo; }

float fo_lms_cc_impl::get_mu_sfo() const { return d_mu_sfo; }

void fo_lms_cc_impl::set_mu_sfo(float mu_sfo) { d_mu_sfo = mu_sfo; }

gr_complex fo_lms_cc_impl::error(const gr_complex& desired, const gr_complex& out) {
    return desired - out;
}

void fo_lms_cc_impl::update_tap(gr_complex& tap, const gr_complex& in) {
    tap += d_mu_cir * in * std::exp(gr_complex(0, 1) * d_p) * std::conj(d_error);
}

// The MMSE interpolator is implemented with 8 taps and, as such, we require that at least 4 samples
// always exist past the currently handled sample.
#define INTERPOLATOR_PADDING 4

int fo_lms_cc_impl::general_work(int noutput_items,
                                 gr_vector_int& ninput_items,
                                 gr_vector_const_void_star& input_items,
                                 gr_vector_void_star& output_items) {
    auto desired = static_cast<const input_type*>(input_items[0]) + d_taps.size() - 1;
    auto input = static_cast<const input_type*>(input_items[1]);
    auto out = static_cast<output_type*>(output_items[0]);
    auto error_out = static_cast<output_type*>(output_items[1]);
    auto cfo_out = static_cast<float*>(output_items[2]);
    auto sfo_out = static_cast<float*>(output_items[3]);

    size_t l = d_taps.size();
    int m = 0;
    int i = 0;

    for (; i < noutput_items - INTERPOLATOR_PADDING; i++) {

        // Increment the phase and time counters.
        d_p += d_cfo;
        if (d_p >= 2 * M_PI) {
            d_p -= 2 * M_PI;
        }
        d_t += (1.0 / d_samp_rate) * d_sfo;
        if (d_t > 1.0) {
            d_t -= 1.0;
            m += 1;
        } else if (d_t < 0.0) {
            d_t += 1.0;
            m -= 1;
        }

        // Perform arbitrary sampling rate conversion.
        for (int k = 0; k < l; k++) {
            // Interpolator gives [0,1] between 4th and 5th sample relying on 8 taps.
            d_y_0[k] = d_interpolator.interpolate(&input[i + k + m - 3], d_t);

            d_y_1[k] = d_interpolator.interpolate(&input[i + k + m - 3 + 1], d_t);
            d_y_2[k] = d_interpolator.interpolate(&input[i + k + m - 3 - 1], d_t);
        }

        // Calculate the output signal y(n) of the adaptive filter.
        volk_32fc_x2_conjugate_dot_prod_32fc(&out[i], &d_y_0[0], &d_taps[0], l);
        out[i] = out[i] * std::exp(gr_complex(0, 1) * d_p);

        // Calculate the error signal e(n) by using: e(n) = d(n) - y(n).
        d_error = error(desired[i], out[i]);

        if (error_out != nullptr) {
            error_out[i] = d_error;
        }

        // Update the filter coefficients.
        for (int k = 0; k < l; k++) {
            // Update tap locally from error.
            update_tap(d_taps[k], d_y_0[k]);
        }

        // Copy previous frequency offset estimates to the respective outputs.
        if (cfo_out != nullptr) {
            cfo_out[i] = d_cfo;
        }
        if (sfo_out != nullptr) {
            sfo_out[i] = d_sfo;
        }

        // Update the carrier frequency offset estimate.
        d_cfo -= d_mu_cfo * std::imag(out[i] * std::conj(d_error));

        // Update the sampling frequency offset estimate.
        gr_complex temp_1, temp_2;
        volk_32fc_x2_conjugate_dot_prod_32fc(&temp_1, &d_y_1[0], &d_taps[0], l);
        volk_32fc_x2_conjugate_dot_prod_32fc(&temp_2, &d_y_2[0], &d_taps[0], l);
        auto d = (temp_1 - temp_2) / gr_complex(2.0, 0);
        d_sfo += d_mu_sfo * std::real(d * std::exp(gr_complex(0, 1) * d_p) * std::conj(d_error));
    }

    // Tell runtime system how many input items we consumed.
    consume(0, i);
    consume(1, i + m);

    // Tell runtime system how many output items we produced.
    return noutput_items - INTERPOLATOR_PADDING;
}

} /* namespace adapt */
} /* namespace gr */
