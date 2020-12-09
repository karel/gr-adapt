/* -*- c++ -*- */
/*
 * Copyright 2018 gr-adapt author.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "iqrd_rls_filter_cc_impl.h"
#include <gnuradio/io_signature.h>
#include <volk/volk.h>

namespace gr {
namespace adapt {

using namespace filter::kernel;

iqrd_rls_filter_cc::sptr iqrd_rls_filter_cc::make(int num_taps,
                                                  float delta,
                                                  float lambda,
                                                  unsigned skip,
                                                  unsigned decimation,
                                                  bool adapt,
                                                  bool reset) {
    return gnuradio::get_initial_sptr(
        new iqrd_rls_filter_cc_impl(num_taps, delta, lambda, skip, decimation, adapt, reset));
}

/*
 * The private constructor
 */
iqrd_rls_filter_cc_impl::iqrd_rls_filter_cc_impl(int num_taps,
                                                 float delta,
                                                 float lambda,
                                                 unsigned skip,
                                                 unsigned decimation,
                                                 bool adapt,
                                                 bool reset)
    : gr::sync_decimator(
          "iqrd_rls_filter_cc",
          gr::io_signature::make(2, 2, sizeof(gr_complex)),
          gr::io_signature::makev(1,
                                  3,
                                  std::vector<int>{sizeof(gr_complex),
                                                   sizeof(gr_complex),
                                                   num_taps * int(sizeof(gr_complex))}),
          decimation),
      fir_filter_ccc(std::vector<gr_complex>(num_taps, gr_complex(0, 0))),
      d_updated(false), d_skip(skip), d_i(0), d_adapt(adapt), d_reset(false) {
    set_delta(delta);
    set_lambda(lambda);

    const int alignment_multiple = volk_get_alignment() / sizeof(gr_complex);
    set_alignment(std::max(1, alignment_multiple));

    set_history(num_taps);

#ifdef ARMADILLO_FOUND
    d_taps.zeros(num_taps);
    d_new_taps.zeros(num_taps);
#else
    d_new_taps = std::vector<gr_complex>(num_taps, gr_complex(0, 0));
#endif // ARMADILLO_FOUND
    init_internals();
}

void iqrd_rls_filter_cc_impl::init_internals() {
    // Soft-start.
#ifdef ARMADILLO_FOUND
    arma::cx_fmat J = d_delta * arma::eye<arma::cx_fmat>(d_taps.size(), d_taps.size());
    d_U = arma::fliplr(J);
    d_dq2.zeros(d_taps.size());
#else
    d_U = std::vector<std::vector<gr_complex>>(
        d_taps.size(), std::vector<gr_complex>(d_taps.size(), gr_complex(0, 0)));
    d_dq2 = std::vector<gr_complex>(d_taps.size(), gr_complex(0.0));
    for (int i = 0; i < d_taps.size(); i++) {
        d_U[d_taps.size() - 1 - i][i] = d_delta;
    }
#endif // ARMADILLO_FOUND
}

void iqrd_rls_filter_cc_impl::set_taps(const std::vector<gr_complex>& new_taps) {
#ifdef ARMADILLO_FOUND
    d_new_taps = arma::cx_fvec(new_taps);
#else
    d_new_taps = new_taps;
#endif // ARMADILLO_FOUND
    d_updated = true;
}

const std::vector<gr_complex>& iqrd_rls_filter_cc_impl::get_taps() {
#ifdef ARMADILLO_FOUND
    fir_filter_ccc::d_taps = arma::conv_to<std::vector<gr_complex>>::from(d_taps);
#endif // ARMADILLO_FOUND
    return fir_filter_ccc::d_taps;
}

float iqrd_rls_filter_cc_impl::get_delta() const { return d_delta; }

void iqrd_rls_filter_cc_impl::set_delta(float delta) {
    if (delta <= 0.0f || delta > 300.0f) {
        throw std::out_of_range(
            "iqrd_rls_filter_cc_impl::set_delta: Regularization factor must be in range (0, 300]");
    } else {
        d_delta = delta;
    }
}

float iqrd_rls_filter_cc_impl::get_lambda() const { return d_lambda; }

void iqrd_rls_filter_cc_impl::set_lambda(float lambda) {
    if (lambda <= 0.0f || lambda > 1.0f) {
        throw std::out_of_range(
            "iqrd_rls_filter_cc_impl::set_lambda: Forgetting factor must be in range (0, 1]");
    } else {
        d_lambda = lambda;
    }
}

unsigned iqrd_rls_filter_cc_impl::get_skip() const { return d_skip; }

void iqrd_rls_filter_cc_impl::set_skip(unsigned skip) { d_skip = skip; }

bool iqrd_rls_filter_cc_impl::get_adapt() const { return d_adapt; }

void iqrd_rls_filter_cc_impl::set_adapt(bool adapt) { d_adapt = adapt; }

bool iqrd_rls_filter_cc_impl::get_reset() const { return d_reset; }

void iqrd_rls_filter_cc_impl::set_reset(bool reset) {
    d_reset = reset;
    if (d_reset) {
        set_taps(std::vector<gr_complex>(d_taps.size(), gr_complex(0, 0)));
    }
}

gr_complex iqrd_rls_filter_cc_impl::error(const gr_complex& desired, const gr_complex& out) {
    return desired - out;
}

int iqrd_rls_filter_cc_impl::work(int noutput_items,
                                  gr_vector_const_void_star& input_items,
                                  gr_vector_void_star& output_items) {
    const auto* desired = (const gr_complex*)input_items[0] + d_taps.size() - 1;
    const auto* input = (const gr_complex*)input_items[1];
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
        init_internals();
        d_updated = false;
        return 0; // history requirements may have changed.
    }

    int k = 0;
    size_t N = d_taps.size();
    gr_complex gamma;
#ifdef ARMADILLO_FOUND
    arma::cx_fvec u;
    arma::cx_fvec input_arma((gr_complex*)input, noutput_items * decimation() + N - 1, false, true);
#endif // ARMADILLO_FOUND
    for (int i = 0; i < noutput_items; i++) {
        if (d_adapt && (!d_skip || d_i >= (d_skip + 1))) {
#ifdef ARMADILLO_FOUND
            // Obtaining a(k)
            arma::cx_fmat akaux(N, 1, arma::fill::zeros);
            arma::cx_fmat xaux =
                (1 / std::sqrt(d_lambda)) * input_arma.subvec(k, arma::size(d_taps));

            for (size_t n = 0; n < N; n++) {
                for (size_t m = 0; m < N - n; m++) {
                    akaux(n) = akaux(n) + (d_U(n, m) * xaux(m));
                }
            }

            arma::cx_fmat a = akaux;

            // Obtaining Q(k) and gamma(k)
            gr_complex igamma(1, 1);
            arma::cx_fmat c(N, 1, arma::fill::zeros);
            arma::cx_fmat s(N, 1, arma::fill::zeros);

            for (size_t n = 0; n < N; n++) {
                gr_complex aux1 = std::hypot(abs(igamma), abs(akaux(N - 1 - n)));
                c(n) = abs(igamma) / aux1;
                s(n) = (akaux(N - 1 - n) / igamma) * c(n);
                // igamma = aux1;
                igamma = c(n) * igamma + std::conj(s(n)) * akaux(N - 1 - n);
            }

            gamma = gr_complex(1, 1) / igamma;

            // Obtaining u(k) and updating U(k)
            arma::cx_fvec uHaux(N, arma::fill::zeros);
            arma::cx_fmat UmHaux = (1 / std::sqrt(d_lambda)) * d_U;

            for (size_t n = 0; n < N; n++) {
                for (size_t m = 0; m <= n; m++) {
                    gr_complex aux2 = uHaux(m);
                    uHaux(m) = c(n) * aux2 - std::conj(s(n)) * UmHaux(N - 1 - n, m);
                    UmHaux(N - 1 - n, m) = s(n) * aux2 + c(n) * UmHaux(N - 1 - n, m);
                }
            }

            u = arma::conj(uHaux);
            d_U = UmHaux;
#else
            // TODO: Implement IQRD-RLS without Armadillo.
#endif // ARMADILLO_FOUND
        }

// Calculate the output signal y(n) of the adaptive filter.
#ifdef ARMADILLO_FOUND
        out[i] = std::conj(arma::dot(input_arma.subvec(k, arma::size(d_taps)).t(), d_taps));
#else
#ifdef ALIGNED_FIR_FILTER
        out[i] = std::conj(filter(&input[k]));
#else
        volk_32fc_x2_dot_prod_32fc(&out[i], &input[k], &d_taps[0], N);
#endif // ALIGNED_FIR_FILTER
#endif // ARMADILLO_FOUND

        // Obtaining e(k)
        gr_complex e = error(desired[k], out[i]);
        if (error_out != nullptr) {
            error_out[i] = e;
        }

        if (d_adapt && (!d_skip || d_i >= (d_skip + 1))) {
            d_i = 0;
#ifdef ARMADILLO_FOUND
            d_taps = d_taps - (gamma * std::conj(e) * u);
#else
            // TODO: Implement taps calculation without Armadillo.
#ifdef ALIGNED_FIR_FILTER
            for (int n = 0; n < N; n++) {
                // Update aligned taps in filter object.
                fir_filter_ccc::update_tap(d_taps[n], n);
            }
#endif // ALIGNED_FIR_FILTER
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

        k += decimation();
    }

    // Tell runtime system how many output items we produced.
    return noutput_items;
}

} /* namespace adapt */
} /* namespace gr */
