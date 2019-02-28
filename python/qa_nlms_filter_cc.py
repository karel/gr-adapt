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

class qa_nlms_filter_cc (gr_unittest.TestCase):
    n_taps = 11
    mu = 1
    skip = 0
    decimation = 1
    adapt = True
    bypass = False
    reset = False
    n_samples = 1024*2
    n_ensemble = 500
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
        nlms_filter_cc = []
        y_sink_nlms = []
        e_sink_nlms = []
        for j, _ in enumerate(W):
            nlms_filter_cc.append(adapt.nlms_filter_cc(True, self.n_taps, self.mu, self.skip, self.decimation, self.adapt, self.bypass, self.reset))
            y_sink_nlms.append(blocks.vector_sink_c())
            e_sink_nlms.append(blocks.vector_sink_c())
            self.tb.connect((nlms_filter_cc[j], 0), y_sink_nlms[j])
            self.tb.connect((nlms_filter_cc[j], 1), e_sink_nlms[j])

        for i in range(0, self.n_ensemble):
            # Set taps to zero
            for j in range(0, len(W)):
                nlms_filter_cc[j].set_taps(np.zeros(self.n_taps, dtype=np.complex128))

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

                self.tb.connect(d_source, (nlms_filter_cc[j], 0))
                self.tb.connect(u_source[j], (nlms_filter_cc[j], 1))

            self.tb.run()

            for j in range(0, len(W)):
                self.tb.disconnect(d_source, (nlms_filter_cc[j], 0))
                self.tb.disconnect(u_source[j], (nlms_filter_cc[j], 1))

        e_data_nlms = []
        for j in range(0, len(W)):
            e_data_nlms.append((np.abs(e_sink_nlms[j].data()) ** 2).reshape((self.n_ensemble, self.n_samples-h[j].size*1)))
            e_data_nlms[j] = np.mean(e_data_nlms[j], axis=0)

        if self.plot_enabled:
            plt.figure(figsize=(5, 4))
            plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k'])))
            plt.title("Learning curves for NLMS algorithm")
            plt.xlabel("Number of iterations")
            plt.ylabel("Ensemble-averaged square error")
            for j in range(0, len(W)):
                plt.semilogy(e_data_nlms[j], label="W={}".format(W[j]))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.show()

    def test_002_t (self):
        '''
        Simulate convergence for different algorithm parameters.
        :return:
        '''
        # Parameters
        mu = [0.9, 1, 1.1]
        n_samples = 1024*2

        # Channel model
        W = [3.1]
        h = np.zeros((len(W), 4))
        for j, Wj in enumerate(W):
            h[j][1:4] = 0.5 * (1 + np.cos(2*np.pi/Wj * np.linspace(-1,1,3)))

        # Filters
        nlms_filter_cc = []
        y_sink_nlms = []
        e_sink_nlms = []
        for j in range(0, len(mu)):
            nlms_filter_cc.append(adapt.nlms_filter_cc(True, self.n_taps, mu[j], self.skip, self.decimation, self.adapt, self.bypass, self.reset))
            y_sink_nlms.append(blocks.vector_sink_c())
            e_sink_nlms.append(blocks.vector_sink_c())
            self.tb.connect((nlms_filter_cc[j], 0), y_sink_nlms[j])
            self.tb.connect((nlms_filter_cc[j], 1), e_sink_nlms[j])

        for i in range(0, self.n_ensemble):
            # Set taps to zero
            for j in range(0, len(mu)):
                nlms_filter_cc[j].set_taps(np.zeros(self.n_taps, dtype=np.float))

            # Some useful signal to be transmitted
            np.random.seed(i)
            d = np.zeros(self.n_samples, dtype=np.complex128)
            d.real = np.random.randint(2, size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence
            d.imag = np.random.randint(2, size=self.n_samples)*2-1
            d_source = blocks.vector_source_c(d.tolist())

            u = []
            u_source = []
            for j in range(0, len(mu)):
                u.append(np.convolve(d, h[0], mode='valid')) # Distorted signal
                u[j] += np.random.normal(0, np.sqrt(0.001), u[j].size) # Received signal

                u_source.append(blocks.vector_source_c(u[j].tolist()))

                self.tb.connect(d_source, (nlms_filter_cc[j], 0))
                self.tb.connect(u_source[j], (nlms_filter_cc[j], 1))

            self.tb.run()

            for j in range(0, len(mu)):
                self.tb.disconnect(d_source, (nlms_filter_cc[j], 0))
                self.tb.disconnect(u_source[j], (nlms_filter_cc[j], 1))

        e_data_nlms = []
        for j in range(0, len(mu)):
            e_data_nlms.append((np.abs(e_sink_nlms[j].data()) ** 2).reshape((self.n_ensemble, n_samples-h[0].size*1)))
            e_data_nlms[j] = np.mean(e_data_nlms[j], axis=0)

        if self.plot_enabled:
            plt.figure(figsize=(5, 4))
            plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k'])))
            plt.title("Learning curves for NLMS algorithm")
            plt.xlabel("Number of iterations")
            plt.ylabel("Ensemble-averaged square error")
            for j in range(0, len(mu)):
                plt.semilogy(e_data_nlms[j], label=u"μ={}".format(mu[j]))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.show()

    def test_003_t (self):
        '''
        Measure performance.
        :return:
        '''
        # Overwrite variables
        self.n_samples = 1024*1000

        # Create signal
        d = np.zeros(self.n_samples, dtype=np.complex128)
        d.real = np.random.randint(2, size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence
        d.imag = np.random.randint(2, size=self.n_samples)*2-1

        # Channel model
        W = 2.9
        h = np.zeros(4)
        h[1:4] = 0.5 * (1 + np.cos(2*np.pi/W * np.linspace(-1,1,3)))

        u = np.convolve(d, h, mode='valid') # Distorted signal
        u += np.random.normal(0, np.sqrt(0.001), u.size) # Received signal

        d_source = blocks.vector_source_c(d.tolist())
        u_source = blocks.vector_source_c(u.tolist())
        nlms_filter_cc = adapt.nlms_filter_cc(True, self.n_taps, self.mu, self.skip, self.decimation, self.adapt, self.bypass, self.reset)
        y_sink = blocks.vector_sink_c()
        e_sink = blocks.vector_sink_c()
        self.tb.connect(d_source, (nlms_filter_cc, 0))
        self.tb.connect(u_source, (nlms_filter_cc, 1))
        self.tb.connect((nlms_filter_cc, 0), y_sink)
        self.tb.connect((nlms_filter_cc, 1), e_sink)
        self.tb.run()
        throughput_avg = nlms_filter_cc.pc_throughput_avg()
        print("pc_throughput_avg: {0:.3f} MSPS".format(throughput_avg / 1e6))

    def test_004_t (self):
        '''
        Check if decimation works.
        :return:
        '''
        # Overwrite variables
        self.n_samples = 1024
        self.decimation = 2

        # Create signal
        Fs = 2e6
        Fc = 100e3

        phaseAcc = 0
        phaseInc = 2 * np.pi * Fc / Fs
        phaseAccNext = phaseAcc + self.n_samples * phaseInc
        d = 1 * np.exp(np.linspace(phaseAcc, phaseAccNext, self.n_samples))

        # Channel model
        W = 3.1
        h = np.zeros(4)
        h[1:4] = 0.5 * (1 + np.cos(2*np.pi/W * np.linspace(-1,1,3)))

        u = np.convolve(d, h, mode='valid') # Distorted signal
        u += np.random.normal(0, np.sqrt(0.001), u.size) # Received signal

        d_source = blocks.vector_source_c(d.tolist())
        u_source = blocks.vector_source_c(u.tolist())
        nlms_filter_cc = adapt.nlms_filter_cc(True, self.n_taps, self.mu, self.skip, self.decimation, self.adapt, self.bypass, self.reset)
        y_sink = blocks.vector_sink_c()
        e_sink = blocks.vector_sink_c()
        self.tb.connect(d_source, (nlms_filter_cc, 0))
        self.tb.connect(u_source, (nlms_filter_cc, 1))
        self.tb.connect((nlms_filter_cc, 0), y_sink)
        self.tb.connect((nlms_filter_cc, 1), e_sink)
        self.tb.run()

        y_data = np.array(y_sink.data())
        e_data = np.abs(e_sink.data())

        plt.figure(figsize=(15, 9))
        plt.subplot(211)
        plt.title("Adaptation")
        plt.xlabel("samples - k")
        plt.plot(np.delete(d, np.arange(0, d.size, 2)), "b", label="d - reference")
        plt.plot(np.delete(u, np.arange(0, u.size, 2)), "r", label="u - input")
        plt.plot(y_data, "g", label="y - output")
        plt.legend()
        plt.subplot(212)
        plt.title("Filter error")
        plt.xlabel("samples - k")
        plt.plot(10 * np.log10(e_data ** 2), "r", label="e - error [dB]")
        plt.legend()
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    gr_unittest.run(qa_nlms_filter_cc, "qa_nlms_filter_cc.xml")
