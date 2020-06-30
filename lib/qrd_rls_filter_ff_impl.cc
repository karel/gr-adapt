/* -*- c++ -*- */
/*
 * Copyright 2018 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qrd_rls_filter_ff_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>

namespace gr {
namespace adapt {

using namespace filter::kernel;

qrd_rls_filter_ff::sptr qrd_rls_filter_ff::make(int num_taps,
                                                float delta,
                                                float lambda,
                                                unsigned skip,
                                                unsigned decimation,
                                                bool adapt,
                                                bool reset) {
    return gnuradio::get_initial_sptr(
        new qrd_rls_filter_ff_impl(num_taps, delta, lambda, skip, decimation, adapt, reset));
}

/*
 * The private constructor
 */
qrd_rls_filter_ff_impl::qrd_rls_filter_ff_impl(int num_taps,
                                               float delta,
                                               float lambda,
                                               unsigned skip,
                                               unsigned decimation,
                                               bool adapt,
                                               bool reset)
    : gr::sync_decimator(
          "qrd_rls_filter_ff",
          gr::io_signature::make(2, 2, sizeof(float)),
          gr::io_signature::makev(
              1, 3, std::vector<int>{sizeof(float), sizeof(float), num_taps * int(sizeof(float))}),
          decimation),
      fir_filter_fff(decimation, std::vector<float>(num_taps, 0.0)), d_updated(false), d_skip(skip),
      d_i(0), d_adapt(adapt), d_reset(false) {
    set_delta(delta);
    set_lambda(lambda);

    const int alignment_multiple = volk_get_alignment() / sizeof(float);
    set_alignment(std::max(1, alignment_multiple));

    set_history(num_taps);

#ifdef ARMADILLO_FOUND
    d_taps.zeros(num_taps);
    d_new_taps.zeros(num_taps);
#else
    d_new_taps = std::vector<float>(num_taps, 0.0);
#endif // ARMADILLO_FOUND
    init_internals();
}

void qrd_rls_filter_ff_impl::init_internals() {
    // Soft-start.
#ifdef ARMADILLO_FOUND
    arma::fmat J = d_delta * arma::eye<arma::fmat>(d_taps.size(), d_taps.size());
    d_U = arma::fliplr(J);
    d_dq2.zeros(d_taps.size());
#else
    d_U = std::vector<std::vector<float>>(d_taps.size(), std::vector<float>(d_taps.size(), 0.0));
    d_dq2 = std::vector<float>(d_taps.size(), 0.0);
    for (int i = 0; i < d_taps.size(); i++) {
        d_U[d_taps.size() - 1 - i][i] = d_delta;
    }
#endif // ARMADILLO_FOUND
}

void qrd_rls_filter_ff_impl::set_taps(const std::vector<float>& new_taps) {
#ifdef ARMADILLO_FOUND
    d_new_taps = arma::fvec(new_taps);
#else
    d_new_taps = new_taps;
#endif // ARMADILLO_FOUND
    d_updated = true;
}

const std::vector<float>& qrd_rls_filter_ff_impl::get_taps() {
#ifdef ARMADILLO_FOUND
    fir_filter_fff::d_taps = arma::conv_to<std::vector<float>>::from(d_taps);
#endif // ARMADILLO_FOUND
    return fir_filter_fff::d_taps;
}

float qrd_rls_filter_ff_impl::get_delta() const { return d_delta; }

void qrd_rls_filter_ff_impl::set_delta(float delta) {
    if (delta <= 0.0f || delta > 300.0f) {
        throw std::out_of_range(
            "qrd_rls_filter_ff_impl::set_delta: Regularization factor must be in range (0, 300]");
    } else {
        d_delta = delta;
    }
}

float qrd_rls_filter_ff_impl::get_lambda() const { return d_lambda; }

void qrd_rls_filter_ff_impl::set_lambda(float lambda) {
    if (lambda <= 0.0f || lambda > 1.0f) {
        throw std::out_of_range(
            "qrd_rls_filter_ff_impl::set_lambda: Forgetting factor must be in range (0, 1]");
    } else {
        d_lambda = lambda;
    }
}

unsigned qrd_rls_filter_ff_impl::get_skip() const { return d_skip; }

void qrd_rls_filter_ff_impl::set_skip(unsigned skip) { d_skip = skip; }

bool qrd_rls_filter_ff_impl::get_adapt() const { return d_adapt; }

void qrd_rls_filter_ff_impl::set_adapt(bool adapt) { d_adapt = adapt; }

bool qrd_rls_filter_ff_impl::get_reset() const { return d_reset; }

void qrd_rls_filter_ff_impl::set_reset(bool reset) {
    d_reset = reset;
    if (d_reset) {
        set_taps(std::vector<float>(d_taps.size(), 0.0));
    }
}

float qrd_rls_filter_ff_impl::error(const float& desired, const float& out) {
    return desired - out;
}

int qrd_rls_filter_ff_impl::work(int noutput_items,
                                 gr_vector_const_void_star& input_items,
                                 gr_vector_void_star& output_items) {
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

    int k = 0;
    size_t N = d_taps.size();
#ifdef ARMADILLO_FOUND
    arma::fvec input_arma((float*)input, noutput_items * decimation() + N - 1, false, true);
#endif // ARMADILLO_FOUND
    for (int i = 0; i < noutput_items; i++) {
        if (d_adapt && (!d_skip || d_i >= (d_skip + 1))) {
            d_i = 0;
#ifdef ARMADILLO_FOUND
            arma::fmat xaux = input_arma.subvec(k, arma::size(d_taps)).t();
            arma::fmat Uaux = sqrtf(d_lambda) * d_U;
            float daux = desired[k];
            arma::fmat dq2aux = sqrtf(d_lambda) * d_dq2;
            float gamma = 1;

            for (size_t n = 0; n < N; n++) {
                // Obtaining Q and updating U
                float c = abs(Uaux(N - 1 - n, n)) / hypotf(xaux(n), Uaux(N - 1 - n, n));
                float s = (xaux(n) / Uaux(N - 1 - n, n)) * c;
                Uaux(N - 1 - n, n) = s * xaux(n) + c * Uaux(N - 1 - n, n);
                xaux(n) = 0;
                for (size_t m = n + 1; m < N; m++) {
                    float oldxaux = xaux(m);
                    xaux(m) = c * oldxaux - s * Uaux(N - 1 - n, m);
                    Uaux(N - 1 - n, m) = s * oldxaux + c * Uaux(N - 1 - n, m);
                }

                // Obtaining gamma
                gamma = gamma * c;

                // Obtaining eq1 and updating dq2
                float olddaux = daux;
                daux = c * olddaux - s * dq2aux(N - 1 - n);
                dq2aux(N - 1 - n) = s * olddaux + c * dq2aux(N - 1 - n);

                // Back substitution if necessary
                float summa = 0;
                for (size_t m = 0; m < n; m++) {
                    summa = summa + Uaux(n, N - 1 - m) * d_taps[N - 1 - m];
                }
                d_taps[N - 1 - n] = (dq2aux(n) - summa) / Uaux(n, N - 1 - n);
#ifdef ALIGNED_FIR_FILTER
                // Update aligned taps in filter object.
                fir_filter_fff::update_tap(d_taps[N - 1 - n], N - 1 - n);
#endif // ALIGNED_FIR_FILTER
            }

            d_U = Uaux;
            d_dq2 = dq2aux;
#else
            // TODO: Implement QRD-RLS without Armadillo.
#endif // ARMADILLO_FOUND
        } else {
            d_i++;
        }

        // Calculate the output signal y(n) of the adaptive filter.
#ifdef ARMADILLO_FOUND
        out[i] = arma::dot(input_arma.subvec(k, arma::size(d_taps)), d_taps);
#else
#ifdef ALIGNED_FIR_FILTER
        out[i] = filter(&input[k]);
#else
        volk_32f_x2_dot_prod_32f(&out[i], &input[k], &d_taps[0], N);
#endif // ALIGNED_FIR_FILTER
#endif // ARMADILLO_FOUND

        if (error_out != nullptr) {
            error_out[i] = error(desired[k], out[i]); // daux / gamma
        }

        if (taps_out != nullptr) {
#ifdef ARMADILLO_FOUND
            std::memcpy(
                &taps_out[i * d_taps.size()], d_taps.memptr(), sizeof(float) * d_taps.size());
#else
            std::memcpy(&taps_out[i * d_taps.size()], &(d_taps[0]), sizeof(float) * d_taps.size());
#endif // ARMADILLO_FOUND
        }

        k += decimation();
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace adapt */
} /* namespace gr */
