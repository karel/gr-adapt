#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2020 gr-adapt author.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import numpy as np
from scipy import special
from gnuradio import gr

class volterra_lms_ff(gr.decim_block):
    """
    docstring for block volterra_lms_ff
    """
    def __init__(self, num_taps, mu1, mu2, adapt, bypass, reset):
        gr.decim_block.__init__(self,
            name="volterra_lms_ff",
            in_sig=[np.float32, np.float32],
            out_sig=[np.float32, np.float32], decim=1)
        self.mu1 = mu1
        self.mu2 = mu2
        self.adapt = adapt
        self.bypass = bypass
        self.reset = reset
        self.N = num_taps # memory
        self.P = 2 # order
        self.length = 0
        for p in range(1, self.P+1):
            self.length += int(special.binom(self.N + p - 1, p))
        self.mus = np.concatenate((np.repeat(mu1, self.N), \
            np.repeat(mu2, self.length-self.N)), axis=None)
        self.x = np.zeros(self.length, dtype=np.float32)
        self.w = np.zeros(self.length, dtype=np.float32)
        self.set_history(self.N)
    
    def set_mu1(self, mu1):
        self.mu1 = mu1
        self.mus = np.concatenate((np.repeat(self.mu1, self.N), \
            np.repeat(self.mu2, self.length-self.N)), axis=None)

    def set_mu2(self, mu2):
        self.mu2 = mu2
        self.mus = np.concatenate((np.repeat(self.mu1, self.N), \
            np.repeat(self.mu2, self.length-self.N)), axis=None)
    
    def set_adapt(self, adapt):
        self.adapt = adapt

    def set_bypass(self, bypass):
        self.bypass = bypass

    def set_reset(self, reset):
        if reset:
            self.w = np.zeros(self.length, dtype=np.float32)

    def work(self, input_items, output_items):
        desired = input_items[0][self.N-1:]
        input = input_items[1]
        out = output_items[0]
        error_out = output_items[1]

        for i in range(0, len(output_items[0])):
            self.x[:self.N] = input[i:i+self.N]
            j = self.N
            for k in range(0, self.N):
                for l in range(k, self.N):
                    self.x[j] = self.x[k] * self.x[l]
                    j += 1        
            out[i] = self.x.T @ self.w
            error_out[i] = desired[i] - out[i]
            
            if self.adapt:
                self.w += self.mus * error_out[i] * self.x

        return len(output_items[0])
