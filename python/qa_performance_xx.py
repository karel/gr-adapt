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

import cpuinfo
from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np
import adapt_swig as adapt
from gnuradio import blocks
from gnuradio import gr, gr_unittest

class qa_performance_xx (gr_unittest.TestCase):
    first_input = True
    n_taps = 11
    mu_lms = 0.075
    mu_nlms = 1.0
    delta_rls = 1
    lambda_rls = 0.99
    delta_qrd_rls = 0.5
    lambda_qrd_rls = 0.9
    delta_iqrd_rls = 0.5
    lambda_iqrd_rls = 0.9
    skip = 0
    decimation = 1
    adapt = True
    bypass = False
    reset = False
    n_samples = 1024*1000
    plot_enabled = True

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        # Create float filters
        filters_ff = []
        lms_filter_ff = adapt.lms_filter_ff(self.first_input, self.n_taps, self.mu_lms, self.skip, self.decimation, self.adapt, self.bypass, self.reset)
        filters_ff.append(lms_filter_ff)

        nlms_filter_ff = adapt.nlms_filter_ff(self.first_input, self.n_taps, self.mu_nlms, self.skip, self.decimation, self.adapt, self.bypass, self.reset)
        filters_ff.append(nlms_filter_ff)

        rls_filter_ff = adapt.rls_filter_ff(self.first_input, self.n_taps, self.delta_rls, self.lambda_rls, self.skip, self.decimation, self.adapt, self.reset)
        filters_ff.append(rls_filter_ff)

        qrd_rls_filter_ff = adapt.qrd_rls_filter_ff(self.n_taps, self.delta_qrd_rls, self.lambda_qrd_rls, self.skip, self.decimation, self.adapt, self.reset)
        filters_ff.append(qrd_rls_filter_ff)

        iqrd_rls_filter_ff = adapt.iqrd_rls_filter_ff(self.n_taps, self.delta_iqrd_rls, self.lambda_iqrd_rls, self.skip, self.decimation, self.adapt, self.reset)
        filters_ff.append(iqrd_rls_filter_ff)

        y_sink_f = blocks.vector_sink_f()
        e_sink_f = blocks.vector_sink_f()

        # Create complex filters
        filters_cc = []
        lms_filter_cc = adapt.lms_filter_cc(self.first_input, self.n_taps, self.mu_lms, self.skip, self.decimation, self.adapt, self.bypass, self.reset)
        filters_cc.append(lms_filter_cc)

        nlms_filter_cc = adapt.nlms_filter_cc(self.first_input, self.n_taps, self.mu_nlms, self.skip, self.decimation, self.adapt, self.bypass, self.reset)
        filters_cc.append(nlms_filter_cc)

        rls_filter_cc = adapt.rls_filter_cc(self.first_input, self.n_taps, self.delta_rls, self.lambda_rls, self.skip, self.decimation, self.adapt, self.reset)
        filters_cc.append(rls_filter_cc)

        qrd_rls_filter_cc = adapt.qrd_rls_filter_cc(self.n_taps, self.delta_qrd_rls, self.lambda_qrd_rls, self.skip, self.decimation, self.adapt, self.reset)
        filters_cc.append(qrd_rls_filter_cc)

        iqrd_rls_filter_cc = adapt.iqrd_rls_filter_cc(self.n_taps, self.delta_iqrd_rls, self.lambda_iqrd_rls, self.skip, self.decimation, self.adapt, self.reset)
        filters_cc.append(iqrd_rls_filter_cc)

        y_sink_c = blocks.vector_sink_c()
        e_sink_c = blocks.vector_sink_c()

        results_f = []
        results_c = []
        names = []

        # Channel model
        W = 3.1
        h = np.zeros(4)
        h[1:4] = 0.5 * (1 + np.cos(2*np.pi/W * np.linspace(-1,1,3)))

        # Some useful signal to be transmitted
        np.random.seed(2701)
        d = np.random.randint(2,size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence

        u = np.convolve(d, h, mode='valid') # Distorted signal
        u += np.random.normal(0, np.sqrt(0.001), u.size) # Received signal

        np.random.seed(2701)
        d_c = np.empty(self.n_samples, dtype=np.complex64)
        d_c.real = np.random.randint(2,size=self.n_samples)*2-1 # Random bipolar (-1,1) sequence
        d_c.imag = np.random.randint(2,size=self.n_samples)*2-1

        u_c = np.convolve(d_c, h, mode='valid') # Distorted signal
        u_c += np.random.normal(0, np.sqrt(0.001), u_c.size) # Received signal

        for filter_ff in filters_ff:
            d_source = blocks.vector_source_f(d.tolist())
            u_source = blocks.vector_source_f(u.tolist())

            self.tb.connect(d_source, (filter_ff, 0))
            self.tb.connect(u_source, (filter_ff, 1))

            self.tb.connect((filter_ff, 0), y_sink_f)
            self.tb.connect((filter_ff, 1), e_sink_f)

            self.tb.run()

            results_f.append(filter_ff.pc_throughput_avg() / 1e6)
            names.append(str(filter_ff.__class__.__name__).replace('_ff_sptr',''))
            self.tb.disconnect(filter_ff)

        for filter_cc in filters_cc:
            d_source = blocks.vector_source_c(d_c.tolist())
            u_source = blocks.vector_source_c(u_c.tolist())

            self.tb.connect(d_source, (filter_cc, 0))
            self.tb.connect(u_source, (filter_cc, 1))

            self.tb.connect((filter_cc, 0), y_sink_c)
            self.tb.connect((filter_cc, 1), e_sink_c)

            self.tb.run()

            results_c.append(filter_cc.pc_throughput_avg() / 1e6)
            # names.append(str(filter_cc.__class__.__name__).replace('_sptr',''))
            self.tb.disconnect(filter_cc)

        if self.plot_enabled:
            # plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k'])))
            fig, ax = plt.subplots(figsize=(5,4))

            index = np.arange(len(names))  # the x locations for the groups
            width = 0.3 # the width of the bars

            rects1 = ax.bar(index, results_f, width, label='Float (n={})'.format(self.n_taps))
            rects2 = ax.bar(index+width, results_c, width, label='Complex (n={})'.format(self.n_taps))

            ax.set_ylabel('Throughput (MSPS)')
            ax.set_title('Throughput for different algorithms\n{}'.format(cpuinfo.get_cpu_info()['brand']))
            ax.set_xticks(index + width / 2)
            ax.set_xticklabels(names, rotation=45, ha="right")
            ax.legend()

            xpos = 'center'
            xpos = xpos.lower()  # normalize the case of the parameter
            ha = {'center': 'center', 'right': 'left', 'left': 'right'}
            offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

            for rect in rects1:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                        '{0:.2f}'.format(height), ha=ha[xpos], va='bottom')

            for rect in rects2:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                        '{0:.2f}'.format(height), ha=ha[xpos], va='bottom')

            fig.tight_layout()
            plt.show()


if __name__ == '__main__':
    gr_unittest.run(qa_performance_xx, "qa_performance_xx.xml")
