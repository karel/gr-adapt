#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2019 <+YOU OR YOUR COMPANY+>.
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

class qa_rls_filter_ff (gr_unittest.TestCase):
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
        rls_filter_ff = []
        y_sink_rls = []
        e_sink_rls = []
        for j, _ in enumerate(W):
            rls_filter_ff.append(adapt.rls_filter_ff(True, self.n_taps, self.delta, self._lambda, self.skip, self.decimation, self.adapt, self.reset))
            y_sink_rls.append(blocks.vector_sink_f())
            e_sink_rls.append(blocks.vector_sink_f())
            self.tb.connect((rls_filter_ff[j], 0), y_sink_rls[j])
            self.tb.connect((rls_filter_ff[j], 1), e_sink_rls[j])

        for i in range(0, self.n_ensemble):
            # Set taps to zero
            for j in range(0, len(W)):
                rls_filter_ff[j].set_taps(np.zeros(self.n_taps, dtype=np.float))

            # Some useful signal to be transmitted
            np.random.seed(i)
            d = np.random.randint(2, size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence
            d_source = blocks.vector_source_f(d.tolist())

            u = []
            u_source = []
            for j in range(0, len(W)):
                u.append(np.convolve(d, h[j], mode='valid')) # Distorted signal
                u[j] += np.random.normal(0, np.sqrt(0.001), u[j].size) # Received signal

                u_source.append(blocks.vector_source_f(u[j].tolist()))

                self.tb.connect(d_source, (rls_filter_ff[j], 0))
                self.tb.connect(u_source[j], (rls_filter_ff[j], 1))

            self.tb.run()

            for j in range(0, len(W)):
                self.tb.disconnect(d_source, (rls_filter_ff[j], 0))
                self.tb.disconnect(u_source[j], (rls_filter_ff[j], 1))

        e_data_rls = []
        for j in range(0, len(W)):
            e_data_rls.append((np.array(e_sink_rls[j].data()) ** 2).reshape((self.n_ensemble, self.n_samples-h[j].size*2)))
            e_data_rls[j] = np.mean(e_data_rls[j], axis=0)

        if self.plot_enabled:
            plt.figure(figsize=(5, 4))
            plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k'])))
            plt.title("Learning curves for RLS algorithm")
            plt.xlabel("Number of iterations")
            plt.ylabel("Ensemble-averaged square error")
            for j in range(0, len(W)):
                plt.semilogy(e_data_rls[j], label="W={}".format(W[j]))
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.show()

    # def test_002_t (self):
    #     '''
    #     Simulate convergence for different algorithm parameters.
    #     :return:
    #     '''
    #     # Parameters
    #     _lambda = [0.0075, 0.025, 0.075]
    #     n_samples = 1024*3
    #
    #     # Channel model
    #     W = [3.1]
    #     h = np.zeros((len(W), 4))
    #     for j, Wj in enumerate(W):
    #         h[j][1:4] = 0.5 * (1 + np.cos(2*np.pi/Wj * np.linspace(-1,1,3)))
    #
    #     # Filters
    #     rls_filter_ff = []
    #     y_sink_rls = []
    #     e_sink_rls = []
    #     for j in range(0, len(_lambda)):
    #         rls_filter_ff.append(adapt.rls_filter_ff(True, self.n_taps, self.delta, _lambda[j], self.decimation, self.adapt, self.reset))
    #         y_sink_rls.append(blocks.vector_sink_f())
    #         e_sink_rls.append(blocks.vector_sink_f())
    #         self.tb.connect((rls_filter_ff[j], 0), y_sink_rls[j])
    #         self.tb.connect((rls_filter_ff[j], 1), e_sink_rls[j])
    #
    #     for i in range(0, self.n_ensemble):
    #         # Set taps to zero
    #         for j in range(0, len(_lambda)):
    #             rls_filter_ff[j].set_taps(np.zeros(self.n_taps, dtype=np.float))
    #
    #         # Some useful signal to be transmitted
    #         np.random.seed(i)
    #         d = np.random.randint(2, size=n_samples)*2-1 # Random bipolar (-1,1) sequence
    #         d_source = blocks.vector_source_f(d.tolist())
    #
    #         u = []
    #         u_source = []
    #         for j in range(0, len(_lambda)):
    #             u.append(np.convolve(d, h[0], mode='valid')) # Distorted signal
    #             u[j] += np.random.normal(0, np.sqrt(0.001), u[j].size) # Received signal
    #
    #             u_source.append(blocks.vector_source_f(u[j].tolist()))
    #
    #             self.tb.connect(d_source, (rls_filter_ff[j], 0))
    #             self.tb.connect(u_source[j], (rls_filter_ff[j], 1))
    #
    #         self.tb.run()
    #
    #         for j in range(0, len(_lambda)):
    #             self.tb.disconnect(d_source, (rls_filter_ff[j], 0))
    #             self.tb.disconnect(u_source[j], (rls_filter_ff[j], 1))
    #
    #     e_data_rls = []
    #     for j in range(0, len(_lambda)):
    #         e_data_rls.append((np.array(e_sink_rls[j].data()) ** 2).reshape((self.n_ensemble, n_samples-h[0].size*2)))
    #         e_data_rls[j] = np.mean(e_data_rls[j], axis=0)
    #
    #     if self.plot_enabled:
    #         plt.figure(figsize=(5, 4))
    #         plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k'])))
    #         plt.title("Learning curves for RLS algorithm")
    #         plt.xlabel("Number of iterations")
    #         plt.ylabel("Ensemble-averaged square error")
    #         for j in range(0, len(_lambda)):
    #             plt.semilogy(e_data_rls[j], label="lambda={}".format(_lambda[j]))
    #         plt.legend()
    #         plt.grid()
    #         plt.tight_layout()
    #         plt.show()

    # def test_003_t (self):
    #     '''
    #     Measure performance.
    #     :return:
    #     '''
    #     # Overwrite variables
    #     self.n_samples = 1024*1000
    #
    #     # Create signal
    #     d = np.random.randint(2, size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence
    #
    #     # Channel model
    #     W = 2.9
    #     h = np.zeros(4)
    #     h[1:4] = 0.5 * (1 + np.cos(2*np.pi/W * np.linspace(-1,1,3)))
    #
    #     u = np.convolve(d, h, mode='valid') # Distorted signal
    #     u += np.random.normal(0, np.sqrt(0.001), u.size) # Received signal
    #
    #     d_source = blocks.vector_source_f(d.tolist())
    #     u_source = blocks.vector_source_f(u.tolist())
    #     rls_filter_ff = adapt.rls_filter_ff(True, self.n_taps, self.delta, self._lambda, self.decimation, self.adapt, self.reset)
    #     y_sink = blocks.vector_sink_f()
    #     e_sink = blocks.vector_sink_f()
    #     self.tb.connect(d_source, (rls_filter_ff, 0))
    #     self.tb.connect(u_source, (rls_filter_ff, 1))
    #     self.tb.connect((rls_filter_ff, 0), y_sink)
    #     self.tb.connect((rls_filter_ff, 1), e_sink)
    #     self.tb.run()
    #     throughput_avg = rls_filter_ff.pc_throughput_avg()
    #     print("pc_throughput_avg: {0:.3f} MSPS".format(throughput_avg / 1e6))


if __name__ == '__main__':
    gr_unittest.run(qa_rls_filter_ff, "qa_rls_filter_ff.xml")
