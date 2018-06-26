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

import matplotlib.pyplot as plt
import numpy as np
from gnuradio import gr, gr_unittest
from gnuradio import blocks
import adapt_swig as adapt

class qa_qrd_rls_filter_ff (gr_unittest.TestCase):
    n_taps = 11
    delta = 0.5
    _lambda = 0.98
    skip = 0
    decimation = 1
    adapt = True
    reset = False
    n_signal = 1024*100
    plot_enabled = False

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        # Create signal
        Fs = 40e6
        Fc = Fs/32

        phaseAcc = 0
        phaseInc = 2 * np.pi * Fc / Fs
        phaseAccNext = phaseAcc + self.n_signal * phaseInc

        d = 1 * np.sin(np.linspace(phaseAcc, phaseAccNext, self.n_signal)).astype(np.float32)
        h = np.array((1, 0.5, 0.3, 0.9)).astype(np.float32)
        u = np.convolve(d, h, mode='valid')
        d = d[h.size-1:]

        d_source = blocks.vector_source_f(d.tolist())
        u_source = blocks.vector_source_f(u.tolist())
        qrd_rls_filter_ff = adapt.qrd_rls_filter_ff(self.n_taps, self.delta, self._lambda, self.skip, self.decimation, self.adapt, self.reset)
        y_sink = blocks.vector_sink_f()
        e_sink = blocks.vector_sink_f()
        self.tb.connect(d_source, (qrd_rls_filter_ff, 0))
        self.tb.connect(u_source, (qrd_rls_filter_ff, 1))
        self.tb.connect((qrd_rls_filter_ff, 0), y_sink)
        self.tb.connect((qrd_rls_filter_ff, 1), e_sink)
        self.tb.run()
        throughput_avg = qrd_rls_filter_ff.pc_throughput_avg()
        print("pc_throughput_avg: {0:.3f} MSPS".format(throughput_avg / 1e6))
        y_data = np.array(y_sink.data())
        e_data = np.array(e_sink.data())
        #self.assertComplexTuplesAlmostEqual(expected_result, result_data, 6)

        if self.plot_enabled:
            plt.figure(figsize=(15, 9))
            plt.subplot(211)
            plt.title("Adaptation")
            plt.xlabel("samples - k")
            plt.plot(d, "b", label="d - reference")
            plt.plot(u, "r", label="u - input")
            plt.plot(y_data, "g", label="y - output")
            plt.legend()
            plt.subplot(212);
            plt.title("Filter error")
            plt.xlabel("samples - k")
            plt.plot(10 * np.log10(e_data ** 2), "r", label="e - error [dB]")
            plt.legend()
            plt.tight_layout()
            plt.show()

    def test_003_t (self):
        # Create noise
        mu, sigma = 0, 0.1 # mean and standard deviation
        d = np.random.normal(mu, sigma, self.n_signal)
        h = np.array((1, 0.5, 0.3, 0.9)).astype(np.float32)
        u = np.convolve(d, h, mode='valid')
        d = d[h.size-1:]

        d_source = blocks.vector_source_f(d.tolist())
        u_source = blocks.vector_source_f(u.tolist())
        qrd_rls_filter_ff = adapt.qrd_rls_filter_ff(self.n_taps, self.delta, self._lambda, self.skip, self.decimation, self.adapt, self.reset)
        y_sink = blocks.vector_sink_f()
        e_sink = blocks.vector_sink_f()
        self.tb.connect(d_source, (qrd_rls_filter_ff, 0))
        self.tb.connect(u_source, (qrd_rls_filter_ff, 1))
        self.tb.connect((qrd_rls_filter_ff, 0), y_sink)
        self.tb.connect((qrd_rls_filter_ff, 1), e_sink)
        self.tb.run()
        throughput_avg = qrd_rls_filter_ff.pc_throughput_avg()
        print("pc_throughput_avg: {0:.3f} MSPS".format(throughput_avg / 1e6))
        y_data = np.array(y_sink.data())
        e_data = np.array(e_sink.data())
        #self.assertComplexTuplesAlmostEqual(expected_result, result_data, 6)

        if self.plot_enabled:
            plt.figure(figsize=(15, 9))
            plt.subplot(211)
            plt.title("Adaptation")
            plt.xlabel("samples - k")
            plt.plot(d, "b", label="d - reference")
            plt.plot(u, "r", label="u - input")
            plt.plot(y_data, "g", label="y - output")
            plt.legend()
            plt.subplot(212);
            plt.title("Filter error")
            plt.xlabel("samples - k")
            plt.plot(10 * np.log10(e_data ** 2), "r", label="e - error [dB]")
            plt.legend()
            plt.tight_layout()
            plt.show()


if __name__ == '__main__':
    gr_unittest.run(qa_qrd_rls_filter_ff, "qa_qrd_rls_filter_xx.xml")
