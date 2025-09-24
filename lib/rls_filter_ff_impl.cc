/* -*- c++ -*- */
/*
 * Copyright 2018 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "rls_filter_ff_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>

namespace gr {
namespace adapt {

using namespace filter::kernel;

rls_filter_ff::sptr rls_filter_ff::make(bool first_input,
                                        int num_taps,
                                        float delta,
                                        float lambda,
                                        unsigned skip,
                                        unsigned decimation,
                                        bool adapt,
                                        bool reset)
{
    return gnuradio::get_initial_sptr(new rls_filter_ff_impl(
        first_input, num_taps, delta, lambda, skip, decimation, adapt, reset));
}

/*
 * The private constructor
 */
rls_filter_ff_impl::rls_filter_ff_impl(bool first_input,
                                       int num_taps,
                                       float delta,
                                       float lambda,
                                       unsigned skip,
                                       unsigned decimation,
                                       bool adapt,
                                       bool reset)
    : gr::sync_decimator(
          "rls_filter_ff",
          gr::io_signature::make(2, 2, sizeof(float)),
          gr::io_signature::makev(1,
                                  3,
                                  std::vector<int>{ sizeof(float),
                                                    sizeof(float),
                                                    num_taps * int(sizeof(float)) }),
          decimation),
      fir_filter_fff(std::vector<float>(num_taps, 0.0)),
      d_first_input(first_input),
      d_updated(false),
      d_delta(0.5),
      d_lambda(0.9),
      d_skip(skip),
      d_i(0),
      d_adapt(adapt),
      d_reset(false)
{
    set_delta(delta);
    set_lambda(lambda);

    const int alignment_multiple = volk_get_alignment() / sizeof(float);
    set_alignment(std::max(1, alignment_multiple));

    set_history(num_taps);

#ifdef ARMADILLO_FOUND
    d_taps.zeros(num_taps);
    d_new_taps.zeros(num_taps);
#else
    d_new_taps = std::vector<float>(d_taps.size(), 0.0);
#endif // ARMADILLO_FOUND
    init_internals();
}

void rls_filter_ff_impl::init_internals()
{
    // Initialize the correlation matrix.
#ifdef ARMADILLO_FOUND
    d_P = (1.0 / d_delta) * arma::eye<arma::fmat>(d_taps.size(), d_taps.size());
#else
    d_P = std::vector<std::vector<float>>(d_taps.size(),
                                          std::vector<float>(d_taps.size(), 0.0));
    for (int i = 0; i < d_taps.size(); i++) {
        d_P[i][i] = 1.0 / d_delta;
    }
#endif // ARMADILLO_FOUND
}

void rls_filter_ff_impl::set_taps(const std::vector<float>& new_taps)
{
#ifdef ARMADILLO_FOUND
    d_new_taps = arma::fvec(new_taps);
#else
    d_new_taps = new_taps;
#endif // ARMADILLO_FOUND
    d_updated = true;
}

const std::vector<float>& rls_filter_ff_impl::get_taps()
{
#ifdef ARMADILLO_FOUND
    fir_filter_fff::d_taps = arma::conv_to<std::vector<float>>::from(d_taps);
#endif // ARMADILLO_FOUND
    return fir_filter_fff::d_taps;
}

float rls_filter_ff_impl::get_delta() const { return d_delta; }

void rls_filter_ff_impl::set_delta(float delta)
{
    if (delta <= 0.0f || delta > 300.0f) {
        throw std::out_of_range("rls_filter_ff_impl::set_delta: Regularization factor "
                                "must be in range (0, 300]");
    } else {
        d_delta = delta;
    }
}

float rls_filter_ff_impl::get_lambda() const { return d_lambda; }

void rls_filter_ff_impl::set_lambda(float lambda)
{
    if (lambda <= 0.0f || lambda > 1.0f) {
        throw std::out_of_range(
            "rls_filter_ff_impl::set_lambda: Forgetting factor must be in range (0, 1]");
    } else {
        d_lambda = lambda;
    }
}

unsigned rls_filter_ff_impl::get_skip() const { return d_skip; }

void rls_filter_ff_impl::set_skip(unsigned skip) { d_skip = skip; }

bool rls_filter_ff_impl::get_adapt() const { return d_adapt; }

void rls_filter_ff_impl::set_adapt(bool adapt) { d_adapt = adapt; }

bool rls_filter_ff_impl::get_reset() const { return d_reset; }

void rls_filter_ff_impl::set_reset(bool reset)
{
    d_reset = reset;
    if (d_reset) {
        set_taps(std::vector<float>(d_taps.size(), 0.0));
    }
}

float rls_filter_ff_impl::error(const float& desired, const float& out)
{
    return desired - out;
}

void rls_filter_ff_impl::update_tap(float& tap, const float& gain)
{
    tap = tap + gain * d_error;
}

int rls_filter_ff_impl::work(int noutput_items,
                             gr_vector_const_void_star& input_items,
                             gr_vector_void_star& output_items)
{
    const auto* desired = (const float*)input_items[0] + d_taps.size() - 1;
    const auto* input = (const float*)input_items[1];
    auto* out = (float*)output_items[0];
    float* error_out = (float*)output_items[1];
    float* taps_out;
    if (output_items.size() == 2) {
        error_out = (float*)output_items[1];
        taps_out = nullptr;
    } else if (output_items.size() == 3) {
        error_out = (float*)output_items[1];
        taps_out = (float*)output_items[2];
    } else {
        error_out = nullptr;
        taps_out = nullptr;
    }

    if (d_updated) {
        d_taps = d_new_taps;
        set_history(d_taps.size());
        init_internals();
        d_updated = false;
        return 0; // history requirements may have changed.
    }

    int j = 0;
    size_t l = d_taps.size();
#ifdef ARMADILLO_FOUND
    arma::fcolvec input_arma(
        (float*)input, noutput_items * decimation() + l - 1, false, true);
#endif // ARMADILLO_FOUND
    for (int i = 0; i < noutput_items; i++) {
        // Calculate the output signal y(n) of the adaptive filter.
#ifdef ARMADILLO_FOUND
        out[i] = arma::dot(input_arma.subvec(j, arma::size(d_taps)), d_taps.t());
#else
#ifdef ALIGNED_FIR_FILTER
        out[i] = filter(&input[j]);
#else
        volk_32f_x2_dot_prod_32f(&out[i], &input[j], &d_taps[0], l);
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

        // Update the filter coefficients.
        if (d_adapt && (!d_skip || d_i >= (d_skip + 1))) {
            d_i = 0;
#ifdef ARMADILLO_FOUND
            arma::fcolvec pi(1 / d_lambda * d_P *
                             input_arma.subvec(j, arma::size(d_taps)));
            float gamma = 1 + 1 / d_lambda *
                                  arma::as_scalar(
                                      input_arma.subvec(j, arma::size(d_taps)).t() * pi);
            arma::fcolvec k(pi / gamma);
            d_taps = d_taps + (k * d_error);
#else
            std::vector<float> pi(l, 0.0);
            float gamma;
            std::vector<float> k(l, 0.0);
            for (int m = 0; m < l; m++) {
                volk_32f_x2_dot_prod_32f(&pi[m], &d_P[m][0], &input[j], l);
                pi[m] = pi[m] / d_lambda;
            }
            volk_32f_x2_dot_prod_32f(&gamma, &input[j], &pi[0], l);
            gamma = 1 + (1 / d_lambda * gamma);
            volk_32f_s32f_multiply_32f(&k[0], &pi[0], 1 / gamma, l);
            for (int m = 0; m < l; m++) {
                // Update tap locally from error.
                update_tap(d_taps[m], k[m]);
#ifdef ALIGNED_FIR_FILTER
                // Update aligned taps in filter object.
                fir_filter_fff::update_tap(d_taps[m], m);
#endif // ALIGNED_FIR_FILTER
            }
#endif // ARMADILLO_FOUND

#ifdef ARMADILLO_FOUND
            d_P = (d_P - k * input_arma.subvec(j, arma::size(d_taps)).t() * d_P) * 1 /
                  d_lambda;
#else
            std::vector<float> P1(l, 0.0);
            for (int m = 0; m < l; m++) {
                for (int n = 0; n < l; n++) {
                    P1[m] += input[j + n] * d_P[n][m];
                }
            }
            for (int m = 0; m < l; m++) {
                for (int n = 0; n < l; n++) {
                    d_P[m][n] = (d_P[m][n] - (k[m] * P1[n])) / d_lambda;
                }
            }
#endif // ARMADILLO_FOUND
        } else {
            d_i++;
        }

        if (taps_out != nullptr) {
#ifdef ARMADILLO_FOUND
            std::memcpy(&taps_out[i * d_taps.size()],
                        d_taps.memptr(),
                        sizeof(float) * d_taps.size());
#else
            std::memcpy(&taps_out[i * d_taps.size()],
                        &(d_taps[0]),
                        sizeof(float) * d_taps.size());
#endif // ARMADILLO_FOUND
        }

        j += decimation();
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace adapt */
} /* namespace gr */
