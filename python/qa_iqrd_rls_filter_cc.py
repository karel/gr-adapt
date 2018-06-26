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
from gnuradio import gr, gr_unittest
from gnuradio import blocks
import adapt_swig as adapt

class qa_iqrd_rls_filter_cc (gr_unittest.TestCase):
    n_taps = 11
    delta = 0.5
    _lambda = 0.99
    skip = 0
    decimation = 1
    adapt = True
    reset = False
    n_samples = 1024*2
    n_ensemble = 100
    plot_enabled = True

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        '''
        Simulate convergence for different channel models.
        :return:
        '''
        # Channel model
        W = [2.9, 3.1, 3.3, 3.5]
        h = np.zeros((len(W), 4))
        for j, Wj in enumerate(W):
            h[j][1:4] = 0.5 * (1 + np.cos(2*np.pi/Wj * np.linspace(-1,1,3)))

        # Filters
        iqrd_rls_filter_cc = []
        y_sink_iqrd_rls = []
        e_sink_iqrd_rls = []
        for j, _ in enumerate(W):
            iqrd_rls_filter_cc.append(adapt.iqrd_rls_filter_cc(self.n_taps, self.delta, self._lambda, self.skip, self.decimation, self.adapt, self.reset))
            y_sink_iqrd_rls.append(blocks.vector_sink_c())
            e_sink_iqrd_rls.append(blocks.vector_sink_c())
            self.tb.connect((iqrd_rls_filter_cc[j], 0), y_sink_iqrd_rls[j])
            self.tb.connect((iqrd_rls_filter_cc[j], 1), e_sink_iqrd_rls[j])

        for i in range(0, self.n_ensemble):
            # Set taps to zero
            for j in range(0, len(W)):
                iqrd_rls_filter_cc[j].set_taps(np.zeros(self.n_taps, dtype=np.complex128))

            # Some useful signal to be transmitted
            np.random.seed(i)
            d = np.zeros(self.n_samples, dtype=np.complex128)
            d.real = np.random.randint(2, size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence
            d.imag = np.random.randint(2, size=self.n_samples)*2-1
            d_source = blocks.vector_source_c(d.tolist())

            u = []
            u_source = []
            for j in range(0, len(W)):
                u.append(np.convolve(d, h[j], mode='valid')) # Distorted signal
                u[j] += np.random.normal(0, np.sqrt(0.001), u[j].size) # Received signal

                u_source.append(blocks.vector_source_c(u[j].tolist()))

                self.tb.connect(d_source, (iqrd_rls_filter_cc[j], 0))
                self.tb.connect(u_source[j], (iqrd_rls_filter_cc[j], 1))

            self.tb.run()

            for j in range(0, len(W)):
                self.tb.disconnect(d_source, (iqrd_rls_filter_cc[j], 0))
                self.tb.disconnect(u_source[j], (iqrd_rls_filter_cc[j], 1))

        e_data_iqrd_rls = []
        for j in range(0, len(W)):
            e_data_iqrd_rls.append((np.abs(e_sink_iqrd_rls[j].data()) ** 2).reshape((self.n_ensemble, self.n_samples-h[j].size*1)))
            e_data_iqrd_rls[j] = np.mean(e_data_iqrd_rls[j], axis=0)

        if self.plot_enabled:
            plt.figure(figsize=(5, 4))
            plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k'])))
            plt.title("Learning curves for IQRD RLS algorithm")
            plt.xlabel("Number of iterations")
            plt.ylabel("Ensemble-averaged square error")
            for j in range(0, len(W)):
                plt.semilogy(e_data_iqrd_rls[j], label="W={}".format(W[j]))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.show()

    def test_002_t (self):
        # Create signal
        Fs = 40e6
        Fc = Fs/32

        phaseAcc = 0
        phaseInc = 2 * np.pi * Fc / Fs
        phaseAccNext = phaseAcc + self.n_samples * phaseInc

        d = 1 * np.exp(1j * np.linspace(phaseAcc, phaseAccNext, self.n_samples)).astype(np.complex64)
        h = np.array((1.0 + 0.4j, 1.0 + 0.3j, 1.0 + 0.4j, 0.2 + 0.3j, 1.0 + 0.4j, 1.0 + 0.3j)).astype(np.complex64)
        u = np.convolve(d, h, mode='valid')
        d = d[h.size-1:]

        d_source = blocks.vector_source_c(d.tolist())
        u_source = blocks.vector_source_c(u.tolist())
        iqrd_rls_filter_cc = adapt.iqrd_rls_filter_cc(self.n_taps, self.delta, self._lambda, self.skip, self.decimation, self.adapt, self.reset)
        y_sink = blocks.vector_sink_c()
        e_sink = blocks.vector_sink_c()
        self.tb.connect(d_source, (iqrd_rls_filter_cc, 0))
        self.tb.connect(u_source, (iqrd_rls_filter_cc, 1))
        self.tb.connect((iqrd_rls_filter_cc, 0), y_sink)
        self.tb.connect((iqrd_rls_filter_cc, 1), e_sink)
        self.tb.run()
        throughput_avg = iqrd_rls_filter_cc.pc_throughput_avg()
        print("pc_throughput_avg: {0:.3f} MSPS".format(throughput_avg / 1e6))
        y_data = np.array(y_sink.data())
        e_data = np.array(e_sink.data())

        if self.plot_enabled:
            plt.figure(figsize=(15, 9))
            plt.subplot(211)
            plt.title("Adaptation")
            plt.xlabel("samples - k")
            plt.plot(d.real, "b", label="d - reference")
            plt.plot(d.imag, "b", label="d - reference")
            plt.plot(u.real, "r", label="u - input")
            plt.plot(u.imag, "r", label="u - input")
            plt.plot(y_data.real, "g", label="y - output")
            plt.plot(y_data.imag, "g", label="y - output")
            plt.legend()
            plt.subplot(212);
            plt.title("Filter error")
            plt.xlabel("samples")
            plt.plot(10 * np.log10(e_data.real ** 2), "r", label="e - error [dB]")
            plt.plot(10 * np.log10(e_data.imag ** 2), "b", label="e - error [dB]")
            plt.legend()
            plt.tight_layout()
            plt.show()


if __name__ == '__main__':
    gr_unittest.run(qa_iqrd_rls_filter_cc, "qa_iqrd_rls_filter_cc.xml")
