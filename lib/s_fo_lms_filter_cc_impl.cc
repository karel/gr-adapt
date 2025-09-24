/* -*- c++ -*- */
/*
 * Copyright 2024 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "s_fo_lms_filter_cc_impl.h"
#include <gnuradio/block.h>
#include <gnuradio/gr_complex.h>
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace gr {
namespace adapt {

using namespace filter::kernel;

using input_type = gr_complex;
using output_type = gr_complex;
s_fo_lms_filter_cc::sptr s_fo_lms_filter_cc::make(double samp_rate,
                                                  int num_taps,
                                                  float mu_cir,
                                                  float mu_cfo,
                                                  float mu_sfo,
                                                  float mu_q,
                                                  bool adapt,
                                                  bool reset)
{
    return gnuradio::make_block_sptr<s_fo_lms_filter_cc_impl>(
        samp_rate, num_taps, mu_cir, mu_cfo, mu_sfo, mu_q, adapt, reset);
}

std::vector<std::vector<float>> transpose(const std::vector<std::vector<float>> X)
{
    std::vector<std::vector<float>> Y(X[0].size(), std::vector<float>(X.size()));
    for (std::vector<float>::size_type i = 0; i < X[0].size(); i++)
        for (std::vector<float>::size_type j = 0; j < X.size(); j++) {
            Y[i][j] = X[j][i];
        }
    return Y;
}

/*
 * The private constructor
 */
s_fo_lms_filter_cc_impl::s_fo_lms_filter_cc_impl(double samp_rate,
                                                 int num_taps,
                                                 float mu_cir,
                                                 float mu_cfo,
                                                 float mu_sfo,
                                                 float mu_q,
                                                 bool adapt,
                                                 bool reset)
    : gr::block("s_fo_lms_filter_cc",
                gr::io_signature::make(2, 2, sizeof(input_type)),
                gr::io_signature::makev(1,
                                        4,
                                        std::vector<int>{ sizeof(output_type),
                                                          sizeof(output_type),
                                                          sizeof(float),
                                                          sizeof(float) })),
      fir_filter_ccc(std::vector<gr_complex>(num_taps, gr_complex(0, 0))),
      d_samp_rate(samp_rate),
      d_M(num_taps),
      d_mu_cir(mu_cir),
      d_mu_cfo(mu_cfo),
      d_mu_sfo(mu_sfo),
      d_mu_q(mu_q),
      d_adapt(adapt),
      d_reset(false)
{
    const int alignment_multiple = volk_get_alignment() / sizeof(input_type);
    set_alignment(std::max(1, alignment_multiple));

    d_y = new gr_complex[d_M]{};
    d_r = new gr_complex[d_M]{};
    std::fill(d_r, d_r + d_M, gr_complex(.0f, .0f));
    d_q = new gr_complex[d_Q]{};
    std::fill(d_q, d_q + d_Q, gr_complex(.0f, .0f));

    for (auto p = 0; p < (d_P + 1); p++) {
        float s = 1.0 / 6.0 * (1.0 / std::pow(d_delta_x, d_P - p));
        std::transform(d_C[p].begin(), d_C[p].end(), d_C[p].begin(), [&s](auto& c) {
            return c * s;
        });
    }

    d_C = transpose(d_C);

    d_S.resize(d_M, std::vector<float>(d_Q, 0.0));
    d_SY.resize(d_Q, std::vector<gr_complex>(d_M, gr_complex(0.0, 0.0)));
    set_history(d_M);
}

s_fo_lms_filter_cc_impl::~s_fo_lms_filter_cc_impl()
{
    delete[] d_y;
    delete[] d_r;
    delete[] d_q;
}

float s_fo_lms_filter_cc_impl::get_mu_cir() const { return d_mu_cir; }

void s_fo_lms_filter_cc_impl::set_mu_cir(float mu_cir) { d_mu_cir = mu_cir; }

float s_fo_lms_filter_cc_impl::get_mu_cfo() const { return d_mu_cfo; }

void s_fo_lms_filter_cc_impl::set_mu_cfo(float mu_cfo) { d_mu_cfo = mu_cfo; }

float s_fo_lms_filter_cc_impl::get_mu_sfo() const { return d_mu_sfo; }

void s_fo_lms_filter_cc_impl::set_mu_sfo(float mu_sfo) { d_mu_sfo = mu_sfo; }

float s_fo_lms_filter_cc_impl::get_mu_q() const { return d_mu_q; }

void s_fo_lms_filter_cc_impl::set_mu_q(float mu_q) { d_mu_q = mu_q; }

bool s_fo_lms_filter_cc_impl::get_adapt() const { return d_adapt; }

void s_fo_lms_filter_cc_impl::set_adapt(bool adapt) { d_adapt = adapt; }

bool s_fo_lms_filter_cc_impl::get_reset() const { return d_reset; }

void s_fo_lms_filter_cc_impl::set_reset(bool reset)
{
    d_reset = reset;
    if (d_reset) {
        d_new_taps = std::vector<gr_complex>(d_taps.size(), gr_complex(0, 0));
        d_updated = true;
    }
}

gr_complex s_fo_lms_filter_cc_impl::error(const gr_complex& desired,
                                          const gr_complex& out)
{
    return desired - out;
}

void s_fo_lms_filter_cc_impl::update_tap(gr_complex& tap, const gr_complex& in)
{
    tap += d_mu_cir * in * std::exp(gr_complex(0, 1) * d_p) * std::conj(d_error);
}

// The MMSE interpolator is implemented with 8 taps and, as such, we require that at least
// 4 samples always exist past the currently handled sample.
#define INTERPOLATOR_PADDING 4

int s_fo_lms_filter_cc_impl::general_work(int noutput_items,
                                          gr_vector_int& ninput_items,
                                          gr_vector_const_void_star& input_items,
                                          gr_vector_void_star& output_items)
{
    auto desired = static_cast<const input_type*>(input_items[0]) + d_taps.size() - 1;
    auto input = static_cast<const input_type*>(input_items[1]);
    auto out = static_cast<output_type*>(output_items[0]);
    auto error_out = static_cast<output_type*>(output_items[1]);
    auto cfo_out = static_cast<float*>(output_items[2]);
    auto sfo_out = static_cast<float*>(output_items[3]);

    if (d_updated) {
        d_cfo = 0;
        d_sfo = 0;
        d_taps = d_new_taps;
        std::fill(d_q, d_q + d_Q, gr_complex(.0f, .0f));
        d_updated = false;
        return 0;
    }

    int m = 0, i = 0, j = 0;
    float u = 0.0;

    for (; i < noutput_items - INTERPOLATOR_PADDING; i++) {
        std::rotate(d_y, d_y + 1, d_y + d_M);
        std::rotate(d_r, d_r + 1, d_r + d_M);

        // Increment the phase and time counters.
        d_p += d_cfo;
        if (d_p >= 2 * M_PIl) {
            d_p -= 2 * M_PIl;
        } else if (d_p <= -2 * M_PIl) {
            d_p += 2 * M_PIl;
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
        d_y[d_M - 1] = d_interpolator.interpolate(&input[i + d_M + m - 3], d_t);

        // Apply the nonlinear model.
        auto a = std::abs(d_y[d_M - 1]) / d_delta_x;
        j = int(std::floor(a)) + 1;
        u = a - j + 1;
        std::vector<float> v{ u * u * u, u * u, u, 1 };

        std::rotate(d_S.begin(), d_S.begin() + 1, d_S.end());
        std::fill(d_S[d_M - 1].begin(), d_S[d_M - 1].end(), 0.0);
        if (j < 1) {
            j += 1;
        } else if (j > 15) {
            j = 15;
        }
        for (int p = 0; p < (d_P + 1); p++)
            volk_32f_x2_dot_prod_32f(
                &d_S[d_M - 1][j - 1 + p], &v[0], &d_C[p][0], (d_P + 1));

        gr_complex temp;
        volk_32fc_32f_dot_prod_32fc(&temp, &d_q[0], &d_S[d_M - 1][0], d_Q);
        d_r[d_M - 1] = d_y[d_M - 1] + d_y[d_M - 1] * temp;

        // Calculate the output signal of the adaptive filter.
        volk_32fc_x2_conjugate_dot_prod_32fc(&out[i], &d_r[0], &d_taps[0], d_M);
        out[i] = out[i] * std::exp(gr_complex(0, 1) * d_p);
        if (i > 0)
            d_out_previous = out[i - 1];

        // Calculate the error signal.
        d_error = error(desired[i], out[i]);

        if (error_out != nullptr) {
            error_out[i] = d_error;
        }

        // Copy previous frequency offset estimates to the respective outputs.
        if (cfo_out != nullptr) {
            cfo_out[i] = d_cfo;
        }
        if (sfo_out != nullptr) {
            sfo_out[i] = d_sfo;
        }

        if (d_adapt) {
            // Update the filter coefficients.
            for (int k = 0; k < d_M; k++) {
                // Update tap locally from error.
                update_tap(d_taps[k], d_r[k]);
            }

            // Update the spline control points.
            for (int k = 0; k < d_M; k++) {
                for (int l = 0; l < d_Q; l++) {
                    d_SY[l][k] = d_S[k][l] * std::conj(d_y[k]);
                }
            }
            for (int l = 0; l < d_Q; l++) {
                gr_complex temporary_accumulator;
                volk_32fc_x2_dot_prod_32fc(
                    &temporary_accumulator, &d_SY[l][0], &d_taps[0], d_M);
                d_q[l] = d_q[l] + d_mu_q * d_error * temporary_accumulator *
                                      std::exp(-gr_complex(0, 1) * d_p);
            }

            // Update the carrier frequency offset estimate.
            d_cfo -= d_mu_cfo * std::imag(out[i] * std::conj(d_error));

            // Update the sampling frequency offset estimate.
            d_sfo += d_mu_sfo * std::real(std::conj(out[i] - d_out_previous) * d_error);
        }
    }

    // Tell runtime system how many input items we consumed.
    consume(0, i);
    consume(1, i + m);

    // Tell runtime system how many output items we produced.
    return noutput_items - INTERPOLATOR_PADDING;
}

} /* namespace adapt */
} /* namespace gr */
