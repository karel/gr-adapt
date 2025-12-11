/* -*- c++ -*- */
/*
 * Copyright 2025 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */


#include "vss_fo_lms_filter_cc_impl.h"
#include <gnuradio/block.h>
#include <gnuradio/gr_complex.h>
#include <gnuradio/io_signature.h>
#include <volk/volk.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <numeric>

namespace gr {
namespace adapt {

using namespace filter::kernel;

using input_type = gr_complex;
using output_type = gr_complex;
vss_fo_lms_filter_cc::sptr
vss_fo_lms_filter_cc::make(double samp_rate,
                           int num_taps,
                           float lambda_e,
                           float lambda_R,
                           float lambda_y,
                           float lambda_epsilon,
                           float lambda_eta,
                           float sigma2_v,
                           const std::vector<float>& step_size_limits,
                           bool adapt,
                           bool reset)
{
    return gnuradio::make_block_sptr<vss_fo_lms_filter_cc_impl>(samp_rate,
                                                                num_taps,
                                                                lambda_e,
                                                                lambda_R,
                                                                lambda_y,
                                                                lambda_epsilon,
                                                                lambda_eta,
                                                                sigma2_v,
                                                                step_size_limits,
                                                                adapt,
                                                                reset);
}


/*
 * The private constructor
 */
vss_fo_lms_filter_cc_impl::vss_fo_lms_filter_cc_impl(
    double samp_rate,
    int num_taps,
    float lambda_e,
    float lambda_R,
    float lambda_y,
    float lambda_epsilon,
    float lambda_eta,
    float sigma2_v,
    const std::vector<float>& step_size_limits,
    bool adapt,
    bool reset)
    : gr::block("vss_fo_lms_filter_cc",
                gr::io_signature::make(2, 2, sizeof(input_type)),
                gr::io_signature::makev(1,
                                        5,
                                        std::vector<int>{ sizeof(output_type),
                                                          sizeof(output_type),
                                                          sizeof(float),
                                                          sizeof(float),
                                                          sizeof(float) })),
      fir_filter_ccc(std::vector<gr_complex>(num_taps, gr_complex(0, 0))),
      d_samp_rate(samp_rate),
      d_mu_cir_min(step_size_limits[0]),
      d_mu_cir_max(step_size_limits[1]),
      d_mu_cfo_min(step_size_limits[2]),
      d_mu_cfo_max(step_size_limits[3]),
      d_mu_sfo_min(step_size_limits[4]),
      d_mu_sfo_max(step_size_limits[5]),
      d_M(num_taps),
      d_lambda_e(lambda_e),
      d_lambda_R(lambda_R),
      d_lambda_y(lambda_y),
      d_lambda_epsilon(lambda_epsilon),
      d_lambda_eta(lambda_eta),
      d_sigma2_v(sigma2_v),
      d_adapt(adapt),
      d_reset(false)
{
    const int alignment_multiple = volk_get_alignment() / sizeof(input_type);
    set_alignment(std::max(1, alignment_multiple));

    d_y = new gr_complex[d_M]{};
    d_Rye = new gr_complex[d_M]{};

    d_mu_cfos = std::vector<float>(d_M, 0.0);
    d_mu_sfos = std::vector<float>(d_M, 0.0);

    set_history(d_M);
}

/*
 * Our virtual destructor.
 */
vss_fo_lms_filter_cc_impl::~vss_fo_lms_filter_cc_impl()
{
    delete[] d_y;
    delete[] d_Rye;
}

float vss_fo_lms_filter_cc_impl::get_lambda_e() const { return d_lambda_e; }

void vss_fo_lms_filter_cc_impl::set_lambda_e(float lambda_e) { d_lambda_e = lambda_e; }

float vss_fo_lms_filter_cc_impl::get_lambda_R() const { return d_lambda_R; }

void vss_fo_lms_filter_cc_impl::set_lambda_R(float lambda_R) { d_lambda_R = lambda_R; }

float vss_fo_lms_filter_cc_impl::get_lambda_y() const { return d_lambda_y; }

void vss_fo_lms_filter_cc_impl::set_lambda_y(float lambda_y) { d_lambda_y = lambda_y; }

float vss_fo_lms_filter_cc_impl::get_lambda_epsilon() const { return d_lambda_epsilon; }

void vss_fo_lms_filter_cc_impl::set_lambda_epsilon(float lambda_epsilon)
{
    d_lambda_epsilon = lambda_epsilon;
}

float vss_fo_lms_filter_cc_impl::get_lambda_eta() const { return d_lambda_eta; }

void vss_fo_lms_filter_cc_impl::set_lambda_eta(float lambda_eta)
{
    d_lambda_eta = lambda_eta;
}

float vss_fo_lms_filter_cc_impl::get_sigma2_v() const { return d_sigma2_v; }

void vss_fo_lms_filter_cc_impl::set_sigma2_v(float sigma2_v) { d_sigma2_v = sigma2_v; }

bool vss_fo_lms_filter_cc_impl::get_adapt() const { return d_adapt; }

void vss_fo_lms_filter_cc_impl::set_adapt(bool adapt) { d_adapt = adapt; }

bool vss_fo_lms_filter_cc_impl::get_reset() const { return d_reset; }

void vss_fo_lms_filter_cc_impl::set_reset(bool reset)
{
    d_reset = reset;
    if (d_reset) {
        d_new_taps = std::vector<gr_complex>(d_taps.size(), gr_complex(0, 0));
        d_updated = true;
    }
}

gr_complex vss_fo_lms_filter_cc_impl::error(const gr_complex& desired,
                                            const gr_complex& out)
{
    return desired - out;
}

void vss_fo_lms_filter_cc_impl::update_tap(gr_complex& tap, const gr_complex& in)
{
    tap += d_mu_cir * in * std::exp(gr_complex(0, 1) * d_p) * std::conj(d_error);
}

// The MMSE interpolator is implemented with 8 taps and, as such, we require that at least
// 4 samples always exist past the currently handled sample.
#define INTERPOLATOR_PADDING 4

int vss_fo_lms_filter_cc_impl::general_work(int noutput_items,
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
    auto noise_out = static_cast<float*>(output_items[4]);

    // Allocate memory for temporary variables
    gr_complex* Rye_forgotten = new gr_complex[d_M]{};
    gr_complex* Rye = new gr_complex[d_M]{};
    float estimate_sigma2_v = 0.0;

    if (d_updated) {
        d_cfo = 0;
        d_sfo = 0;
        d_p = 0.0;
        d_t = 0.0;
        d_sigma2_e = 1.0;
        d_sigma2_y = 1.0;
        std::memset((void*)d_Rye, 0, d_M * sizeof(gr_complex));
        d_taps = d_new_taps;
        d_updated = false;
        return 0;
    }

    int m = 0, i = 0;

    for (; i < noutput_items - INTERPOLATOR_PADDING; i++) {
        std::rotate(d_y, d_y + 1, d_y + d_M);

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

        // Calculate the output signal y(n) of the adaptive filter.
        volk_32fc_x2_conjugate_dot_prod_32fc(&out[i], &d_y[0], &d_taps[0], d_M);
        out[i] = out[i] * std::exp(gr_complex(0, 1) * d_p);
        if (i > 0)
            d_out_previous = out[i - 1];

        // Calculate the error signal e(n) by using: e(n) = d(n) - y(n).
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

        // Update running averages of the signal power levels.
        d_sigma2_e =
            d_lambda_e * d_sigma2_e + (1 - d_lambda_e) * std::pow(std::abs(d_error), 2);
        d_sigma2_y = d_lambda_y * d_sigma2_y +
                     (1 - d_lambda_y) * std::pow(std::abs(d_y[d_M - 1]), 2);

        // Use known measurement noise level or estimate it.
        if (d_sigma2_v)
            estimate_sigma2_v = d_sigma2_v;
        else {
            gr_complex lambda_R = gr_complex(d_lambda_R, 0);
            volk_32fc_s32fc_multiply2_32fc(Rye_forgotten, d_Rye, &lambda_R, d_M);

            gr_complex t =
                (1 - d_lambda_R) * std::exp(gr_complex(0, 1) * d_p) * std::conj(d_error);
            volk_32fc_s32fc_multiply2_32fc(Rye, d_y, &t, d_M);
            volk_32fc_x2_add_32fc(d_Rye, Rye_forgotten, Rye, d_M);

            gr_complex Rye_norm;
            volk_32fc_x2_conjugate_dot_prod_32fc(&Rye_norm, Rye, Rye, d_M);
            estimate_sigma2_v = d_sigma2_e - 1 / d_sigma2_y * std::abs(Rye_norm);

            // Hardcoded lower limit to avoid instability.
            if (estimate_sigma2_v < 1e-12)
                estimate_sigma2_v = 1e-12;
        }

        if (noise_out != nullptr) {
            noise_out[i] = estimate_sigma2_v;
        }

        gr_complex y_norm;
        volk_32fc_x2_conjugate_dot_prod_32fc(&y_norm, d_y, d_y, d_M);
        d_mu_cir =
            (1 - (std::sqrt(estimate_sigma2_v) / (std::sqrt(d_sigma2_e) + 1e-20))) /
            std::abs(y_norm);

        // Limit step size to improve stability.
        if (d_mu_cir > d_mu_cir_max)
            d_mu_cir = d_mu_cir_max;
        else if (d_mu_cir < d_mu_cir_min)
            d_mu_cir = d_mu_cir_min;

        gr_complex h_norm;
        volk_32fc_x2_conjugate_dot_prod_32fc(&h_norm, &d_taps[0], &d_taps[0], d_M);

        float g_cfo = std::imag(out[i] * std::conj(d_error));
        float g_sfo = std::real(std::conj(out[i] - d_out_previous) * d_error);

        d_g_cfo = d_lambda_epsilon * d_g_cfo + (1 - d_lambda_epsilon) * g_cfo;
        d_g_sfo = d_lambda_eta * d_g_sfo + (1 - d_lambda_eta) * g_sfo;

        float d_mu_cfos_mean =
            std::accumulate(d_mu_cfos.begin(), d_mu_cfos.end(), 0.0) / d_mu_cfos.size();
        float b = d_g_cfo * d_mu_cfos_mean;

        float d_mu_sfos_mean =
            std::accumulate(d_mu_sfos.begin(), d_mu_sfos.end(), 0.0) / d_mu_sfos.size();
        float c = d_g_sfo * d_mu_sfos_mean;

        // Calculate variable step sizes for frequency offset updates.
        d_mu_cfo =
            8 *
            std::pow(std::pow(b, 2) * d_mu_cir /
                         (2 * std::pow(std::abs(h_norm), 2) * d_mu_cir *
                              estimate_sigma2_v * std::pow(d_sigma2_y, 2) +
                          std::pow(std::abs(h_norm), 2) * estimate_sigma2_v * d_sigma2_y +
                          1e-10),
                     1.0 / 3.0);

        d_mu_sfo = std::pow(
            std::pow(c, 2) * d_mu_cir /
                (2 * std::pow(std::abs(h_norm), 2) * d_mu_cir * estimate_sigma2_v *
                     std::pow(d_sigma2_y, 2) +
                 std::pow(std::abs(h_norm), 2) * estimate_sigma2_v * d_sigma2_y + 1e-10),
            1.0 / 3.0);

        // Limit step sizes to improve stability.
        if (d_mu_cfo > d_mu_cfo_max)
            d_mu_cfo = d_mu_cfo_max;
        else if (d_mu_cfo < d_mu_cfo_min)
            d_mu_cfo = d_mu_cfo_min;

        if (d_mu_sfo > d_mu_sfo_max)
            d_mu_sfo = d_mu_sfo_max;
        else if (d_mu_sfo < d_mu_sfo_min)
            d_mu_sfo = d_mu_sfo_min;

        std::rotate(d_mu_cfos.begin(), d_mu_cfos.begin() + 1, d_mu_cfos.end());
        std::rotate(d_mu_sfos.begin(), d_mu_sfos.begin() + 1, d_mu_sfos.end());

        d_mu_cfos[d_M - 1] = d_mu_cfo;
        d_mu_sfos[d_M - 1] = d_mu_sfo;

        if (d_adapt) {
            // Update the filter coefficients.
            for (int k = 0; k < d_M; k++) {
                // Update tap locally from error.
                update_tap(d_taps[k], d_y[k]);
            }

            // Update the carrier frequency offset estimate.
            d_cfo -= d_mu_cfo * std::imag(out[i] * std::conj(d_error));

            // Update the sampling frequency offset estimate.
            d_sfo += d_mu_sfo * std::real(std::conj(out[i] - d_out_previous) * d_error);
        }
    }

    //  Free allocated memory.
    delete[] Rye;
    delete[] Rye_forgotten;

    // Tell runtime system how many input items we consumed.
    consume(0, i);
    consume(1, i + m);

    // Tell runtime system how many output items we produced.
    return noutput_items - INTERPOLATOR_PADDING;
}

} /* namespace adapt */
} /* namespace gr */
