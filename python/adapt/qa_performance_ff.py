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

import cpuinfo
from cycler import cycler
import matplotlib.pyplot as plt
import numpy as np
import adapt_swig as adapt
from gnuradio import blocks
from gnuradio import gr, gr_unittest

class qa_performance_ff (gr_unittest.TestCase):
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
        # Create filters
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

        y_sink = blocks.vector_sink_f()
        e_sink = blocks.vector_sink_f()

        results = []
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

        for filter_ff in filters_ff:
            d_source = blocks.vector_source_f(d.tolist())
            u_source = blocks.vector_source_f(u.tolist())

            self.tb.connect(d_source, (filter_ff, 0))
            self.tb.connect(u_source, (filter_ff, 1))

            self.tb.connect((filter_ff, 0), y_sink)
            self.tb.connect((filter_ff, 1), e_sink)

            self.tb.run()

            results.append(filter_ff.pc_throughput_avg() / 1e6)
            names.append(str(filter_ff.__class__.__name__).replace('_sptr',''))
            self.tb.disconnect(filter_ff)

        if self.plot_enabled:
            # plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'k'])))
            fig, ax = plt.subplots(figsize=(5,4))

            index = np.arange(len(names))  # the x locations for the groups
            width = 0.6 # the width of the bars

            rects = ax.bar(index+width/2, results, width, label='Float (n={})'.format(self.n_taps))

            ax.set_ylabel('Throughput (MSPS)')
            ax.set_title('Throughput for different algorithms\n{}'.format(cpuinfo.get_cpu_info()['brand']))
            ax.set_xticks(index + width / 2)
            ax.set_xticklabels(names, rotation=45, ha="right")
            ax.legend()

            xpos = 'center'
            xpos = xpos.lower()  # normalize the case of the parameter
            ha = {'center': 'center', 'right': 'left', 'left': 'right'}
            offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                        '{0:.2f}'.format(height), ha=ha[xpos], va='bottom')

            fig.tight_layout()
            plt.show()
            plt.savefig(fname='qa_performance_ff.png')


if __name__ == '__main__':
    gr_unittest.run(qa_performance_ff, "qa_performance_ff.xml")
