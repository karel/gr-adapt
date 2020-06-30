/* -*- c++ -*- */
/*
 * Copyright 2018 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "nlms_filter_cc_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>

namespace gr {
namespace adapt {

using namespace filter::kernel;

nlms_filter_cc::sptr nlms_filter_cc::make(bool first_input,
                                          int num_taps,
                                          float mu,
                                          unsigned skip,
                                          unsigned decimation,
                                          bool adapt,
                                          bool bypass,
                                          bool reset) {
    return gnuradio::get_initial_sptr(
        new nlms_filter_cc_impl(first_input, num_taps, mu, skip, decimation, adapt, bypass, reset));
}

/*
 * The private constructor
 */
nlms_filter_cc_impl::nlms_filter_cc_impl(bool first_input,
                                         int num_taps,
                                         float mu,
                                         unsigned skip,
                                         unsigned decimation,
                                         bool adapt,
                                         bool bypass,
                                         bool reset)
    : gr::sync_decimator(
          "nlms_filter_cc",
          gr::io_signature::make(2, 3, sizeof(gr_complex)),
          gr::io_signature::makev(1,
                                  3,
                                  std::vector<int>{sizeof(gr_complex),
                                                   sizeof(gr_complex),
                                                   num_taps * int(sizeof(gr_complex))}),
          decimation),
      fir_filter_ccc(decimation, std::vector<gr_complex>(num_taps, gr_complex(0, 0))),
      d_first_input(first_input), d_updated(false),
      d_epsilon(std::numeric_limits<float>::epsilon()), d_skip(skip), d_i(0), d_adapt(adapt),
      d_bypass(bypass), d_reset(false) {
    set_mu(mu);

    const int alignment_multiple = volk_get_alignment() / sizeof(gr_complex);
    set_alignment(std::max(1, alignment_multiple));

    set_history(num_taps);

#ifdef ARMADILLO_FOUND
    d_taps.zeros(num_taps);
    d_new_taps.zeros(num_taps);
#else
    d_new_taps = std::vector<gr_complex>(num_taps, gr_complex(0, 0));
#endif // ARMADILLO_FOUND
}

void nlms_filter_cc_impl::set_taps(const std::vector<gr_complex>& new_taps) {
#ifdef ARMADILLO_FOUND
    d_new_taps = arma::cx_fvec(new_taps);
#else
    d_new_taps = new_taps;
#endif // ARMADILLO_FOUND
    d_updated = true;
}

const std::vector<gr_complex>& nlms_filter_cc_impl::get_taps() {
#ifdef ARMADILLO_FOUND
    fir_filter_ccc::d_taps = arma::conv_to<std::vector<gr_complex>>::from(d_taps);
#endif // ARMADILLO_FOUND
    return fir_filter_ccc::d_taps;
}

float nlms_filter_cc_impl::get_mu() const { return d_mu; }

void nlms_filter_cc_impl::set_mu(float mu) {
    if (mu <= 0.0f || mu >= 2.0f) {
        throw std::out_of_range("nlms_filter_cc_impl::set_mu: Step size must be in range (0, 2)");
    } else {
        d_mu = mu;
    }
}

unsigned nlms_filter_cc_impl::get_skip() const { return d_skip; }

void nlms_filter_cc_impl::set_skip(unsigned skip) { d_skip = skip; }

bool nlms_filter_cc_impl::get_adapt() const { return d_adapt; }

void nlms_filter_cc_impl::set_adapt(bool adapt) { d_adapt = adapt; }

bool nlms_filter_cc_impl::get_bypass() const { return d_bypass; }

void nlms_filter_cc_impl::set_bypass(bool bypass) { d_bypass = bypass; }

bool nlms_filter_cc_impl::get_reset() const { return d_reset; }

void nlms_filter_cc_impl::set_reset(bool reset) {
    d_reset = reset;
    if (d_reset) {
        set_taps(std::vector<gr_complex>(d_taps.size(), gr_complex(0, 0)));
    }
}

gr_complex nlms_filter_cc_impl::error(const gr_complex& desired, const gr_complex& out) {
    return desired - out;
}

void nlms_filter_cc_impl::update_tap(gr_complex& tap, const gr_complex& in) {
    tap += (d_mu / d_power) * conj(in) * d_error;
}

int nlms_filter_cc_impl::work(int noutput_items,
                              gr_vector_const_void_star& input_items,
                              gr_vector_void_star& output_items) {
    const auto* desired = (const gr_complex*)input_items[0] + d_taps.size() - 1;
    const auto* input = (const gr_complex*)input_items[1];
    const gr_complex* filtered_input;
    if (input_items.size() == 3) {
        filtered_input = (gr_complex*)input_items[2];
    } else {
        filtered_input = (gr_complex*)input_items[1];
    }
    auto* out = (gr_complex*)output_items[0];
    gr_complex* error_out;
    gr_complex* taps_out;
    if (output_items.size() == 2) {
        error_out = (gr_complex*)output_items[1];
        taps_out = nullptr;
    } else if (output_items.size() == 3) {
        error_out = (gr_complex*)output_items[1];
        taps_out = (gr_complex*)output_items[2];
    } else {
        error_out = nullptr;
        taps_out = nullptr;
    }

    if (d_updated) {
        d_taps = d_new_taps;
        set_history(d_taps.size());
        d_updated = false;
        return 0; // history requirements may have changed.
    }

    if (d_bypass) {
        std::memcpy(out, input + d_taps.size() - 1, sizeof(gr_complex) * noutput_items);
        if (error_out != nullptr) {
            std::memset(error_out, 0, sizeof(gr_complex) * noutput_items);
        }
        if (taps_out != nullptr) {
            std::memset(taps_out, 0, sizeof(gr_complex) * noutput_items * d_taps.size());
        }
        return noutput_items;
    }

    int j = 0;
    size_t l = d_taps.size();
#ifdef ARMADILLO_FOUND
    gr_complex scale;
    arma::cx_fvec input_arma((gr_complex*)input, noutput_items * decimation() + l - 1, false, true);
    arma::cx_fvec filtered_input_arma(
        (gr_complex*)filtered_input, noutput_items * decimation() + l - 1, false, true);
#endif // ARMADILLO_FOUND
    for (int i = 0; i < noutput_items; i++) {
        // Calculate the output signal y(n) of the adaptive filter.
#ifdef ARMADILLO_FOUND
        out[i] = arma::dot(input_arma.subvec(j, arma::size(d_taps)), d_taps);
#else
#ifdef ALIGNED_FIR_FILTER
        out[i] = filter(&input[j]);
#else
        volk_32fc_x2_dot_prod_32fc(&out[i], &input[j], &d_taps[0], l);
#endif // ALIGNED_FIR_FILTER
#endif // ARMADILLO_FOUND

        // Calculate the error signal e(n) by using: e(n) = d(n) - y(n).
        if (d_first_input == true) { // First input is the reference signal
            d_error = error(desired[j], out[i]);
        } else { // First input is the error signal
            d_error = desired[j];
        }

        if (error_out != nullptr) {
            error_out[i] = d_error;
        }

        if (d_adapt && (!d_skip || d_i >= (d_skip + 1))) {
            d_i = 0;
            // Calculate the power.
#ifdef ARMADILLO_FOUND
            d_power = arma::cdot(filtered_input_arma.subvec(j, arma::size(d_taps)),
                                 filtered_input_arma.subvec(j, arma::size(d_taps)));
#else
            volk_32fc_x2_conjugate_dot_prod_32fc(&d_power, &filtered_input[j], &filtered_input[j], (unsigned)l);
#endif // ARMADILLO_FOUND
            d_power += d_epsilon;

            // Update the filter coefficients.
#ifdef ARMADILLO_FOUND
            scale = d_mu * d_error / d_power;
            d_taps += arma::conj(filtered_input_arma.subvec(j, arma::size(d_taps))) * scale;
#else
            for (int k = 0; k < l; k++) {
                // Update tap locally from error.
                update_tap(d_taps[k], filtered_input[j + k]);
#ifdef ALIGNED_FIR_FILTER
                // Update aligned taps in filter object.
                fir_filter_ccc::update_tap(d_taps[k], k);
#endif // ALIGNED_FIR_FILTER
            }
#endif // ARMADILLO_FOUND
        } else {
            d_i++;
        }

        if (taps_out != nullptr) {
#ifdef ARMADILLO_FOUND
            std::memcpy(
                &taps_out[i * d_taps.size()], d_taps.memptr(), sizeof(gr_complex) * d_taps.size());
#else
            std::memcpy(
                &taps_out[i * d_taps.size()], &(d_taps[0]), sizeof(gr_complex) * d_taps.size());
#endif // ARMADILLO_FOUND
        }

        j += decimation();
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace adapt */
} /* namespace gr */
