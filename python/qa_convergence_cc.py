#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# Copyright 2018 <+YOU OR YOUR COMPANY+>.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np
import adapt_swig as adapt
from gnuradio import blocks
from gnuradio import gr, gr_unittest

class qa_convergence_cc (gr_unittest.TestCase):
    first_input = True
    n_taps = 11
    mu_lms = 0.0075
    mu_nlms = 1.0
    delta_rls = 0.5
    lambda_rls = 0.9
    delta_qrd_rls = 0.5
    lambda_qrd_rls = 0.9
    delta_iqrd_rls = 0.5
    lambda_iqrd_rls = 0.9
    skip = 0
    decimation = 1
    adapt = True
    reset = False
    n_samples = 1024/2
    n_ensemble = 500
    plot_enabled = True

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        # Create filters
        lms_filter_cc = adapt.lms_filter_cc(self.first_input, self.n_taps, self.mu_lms, self.skip, self.decimation, self.adapt, self.reset)
        y_sink_lms = blocks.vector_sink_c()
        e_sink_lms = blocks.vector_sink_c()
        self.tb.connect((lms_filter_cc, 0), y_sink_lms)
        self.tb.connect((lms_filter_cc, 1), e_sink_lms)

        nlms_filter_cc = adapt.nlms_filter_cc(self.first_input, self.n_taps, self.mu_nlms, self.skip, self.decimation, self.adapt, self.reset)
        y_sink_nlms = blocks.vector_sink_c()
        e_sink_nlms = blocks.vector_sink_c()
        self.tb.connect((nlms_filter_cc, 0), y_sink_nlms)
        self.tb.connect((nlms_filter_cc, 1), e_sink_nlms)

        rls_filter_cc = adapt.rls_filter_cc(self.first_input, self.n_taps, self.delta_rls, self.lambda_rls, self.skip, self.decimation, self.adapt, self.reset)
        y_sink_rls = blocks.vector_sink_c()
        e_sink_rls = blocks.vector_sink_c()
        self.tb.connect((rls_filter_cc, 0), y_sink_rls)
        self.tb.connect((rls_filter_cc, 1), e_sink_rls)

        qrd_rls_filter_cc = adapt.qrd_rls_filter_cc(self.n_taps, self.delta_qrd_rls, self.lambda_qrd_rls, self.skip, self.decimation, self.adapt, self.reset)
        y_sink_qrd_rls = blocks.vector_sink_c()
        e_sink_qrd_rls = blocks.vector_sink_c()
        self.tb.connect((qrd_rls_filter_cc, 0), y_sink_qrd_rls)
        self.tb.connect((qrd_rls_filter_cc, 1), e_sink_qrd_rls)

        iqrd_rls_filter_cc = adapt.iqrd_rls_filter_cc(self.n_taps, self.delta_iqrd_rls, self.lambda_iqrd_rls, self.skip, self.decimation, self.adapt, self.reset)
        y_sink_iqrd_rls = blocks.vector_sink_c()
        e_sink_iqrd_rls = blocks.vector_sink_c()
        self.tb.connect((iqrd_rls_filter_cc, 0), y_sink_iqrd_rls)
        self.tb.connect((iqrd_rls_filter_cc, 1), e_sink_iqrd_rls)

        # Channel model
        W = 3.1
        h = np.zeros(4)
        h[1:4] = 0.5 * (1 + np.cos(2*np.pi/W * np.linspace(-1,1,3)))

        for i in range(0, self.n_ensemble):
            # Set taps to zero
            lms_filter_cc.set_taps(np.zeros(self.n_taps, dtype=np.complex128))
            nlms_filter_cc.set_taps(np.zeros(self.n_taps, dtype=np.complex128))
            rls_filter_cc.set_taps(np.zeros(self.n_taps, dtype=np.complex128))
            qrd_rls_filter_cc.set_taps(np.zeros(self.n_taps, dtype=np.complex128))
            iqrd_rls_filter_cc.set_taps(np.zeros(self.n_taps, dtype=np.complex128))

            # Some useful signal to be transmitted
            np.random.seed(i)
            d = np.zeros(self.n_samples, dtype=np.complex128)
            d.real = np.random.randint(2, size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence
            d.imag = np.random.randint(2, size=self.n_samples)*2-1

            u = np.convolve(d, h, mode='valid') # Distorted signal
            u += np.random.normal(0, np.sqrt(0.001), u.size) # Received signal

            d_source = blocks.vector_source_c(d.tolist())
            u_source = blocks.vector_source_c(u.tolist())

            self.tb.connect(d_source, (lms_filter_cc, 0))
            self.tb.connect(u_source, (lms_filter_cc, 1))

            self.tb.connect(d_source, (nlms_filter_cc, 0))
            self.tb.connect(u_source, (nlms_filter_cc, 1))

            self.tb.connect(d_source, (rls_filter_cc, 0))
            self.tb.connect(u_source, (rls_filter_cc, 1))

            self.tb.connect(d_source, (qrd_rls_filter_cc, 0))
            self.tb.connect(u_source, (qrd_rls_filter_cc, 1))

            self.tb.connect(d_source, (iqrd_rls_filter_cc, 0))
            self.tb.connect(u_source, (iqrd_rls_filter_cc, 1))

            self.tb.run()

            self.tb.disconnect(d_source, (lms_filter_cc, 0))
            self.tb.disconnect(u_source, (lms_filter_cc, 1))

            self.tb.disconnect(d_source, (nlms_filter_cc, 0))
            self.tb.disconnect(u_source, (nlms_filter_cc, 1))

            self.tb.disconnect(d_source, (rls_filter_cc, 0))
            self.tb.disconnect(u_source, (rls_filter_cc, 1))

            self.tb.disconnect(d_source, (qrd_rls_filter_cc, 0))
            self.tb.disconnect(u_source, (qrd_rls_filter_cc, 1))

            self.tb.disconnect(d_source, (iqrd_rls_filter_cc, 0))
            self.tb.disconnect(u_source, (iqrd_rls_filter_cc, 1))

        e_data_lms = (np.abs(e_sink_lms.data()) ** 2).reshape((self.n_ensemble, self.n_samples-h.size))
        e_data_nlms = (np.abs(e_sink_nlms.data()) ** 2).reshape((self.n_ensemble, self.n_samples-h.size))
        e_data_rls = (np.abs(e_sink_rls.data()) ** 2).reshape((self.n_ensemble, self.n_samples-h.size))
        e_data_qrd_rls = (np.abs(e_sink_qrd_rls.data()) ** 2).reshape((self.n_ensemble, self.n_samples-h.size))
        e_data_iqrd_rls = (np.abs(e_sink_iqrd_rls.data()) ** 2).reshape((self.n_ensemble, self.n_samples-h.size))

        e_data_lms = np.mean(e_data_lms, axis=0)
        e_data_nlms = np.mean(e_data_nlms, axis=0)
        e_data_rls = np.mean(e_data_rls, axis=0)
        e_data_qrd_rls = np.mean(e_data_qrd_rls, axis=0)
        e_data_iqrd_rls = np.mean(e_data_iqrd_rls, axis=0)

        if self.plot_enabled:
            plt.figure(figsize=(5, 4))
            plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k', 'm'])))
            plt.title("Learning curves for different algorithms")
            plt.xlabel("Number of iterations")
            plt.ylabel("Ensemble-averaged square error")
            plt.semilogy(e_data_lms, label=u"LMS (μ={})".format(self.mu_lms))
            plt.semilogy(e_data_nlms, label=u"NLMS (μ={})".format(self.mu_nlms))
            plt.semilogy(e_data_rls, label=u"RLS (δ={}, λ={})".format(self.delta_rls, self.lambda_rls))
            plt.semilogy(e_data_qrd_rls, label=u"QRD RLS (δ={}, λ={})".format(self.delta_qrd_rls, self.lambda_qrd_rls))
            plt.semilogy(e_data_iqrd_rls, label=u"IQRD RLS (δ={}, λ={})".format(self.delta_iqrd_rls, self.lambda_iqrd_rls))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.show()


if __name__ == '__main__':
    gr_unittest.run(qa_convergence_cc, "qa_convergence_cc.xml")
