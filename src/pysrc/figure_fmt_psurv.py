#!/usr/bin/env python3

##############################################
# name: figure_fmt_psurv.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: probability of mutant bottleneck
# survival as function of degree of escape
##############################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
from scipy.special import lambertw

import plotting_style as ps
import parse_parameters as pp

from analytic_within_host import analytic_p_inocs, neutralize_prob_from_z_poisson, z_poisson

from plotting_functions import add_bounding_subplot

mpl.rcParams['font.size'] = 25

def plot_fmt_psurv(axis = None,
                   bottleneck = None,
                   inoculum_size = None,
                   wt_neut_prob = None,
                   cross_immunities = None,
                   legend = True,
                   lw = 5,
                   plot_drift = True,
                   crop_low = 0.4,
                   crop_high = 0.8,
                   fineness = 50,
                   max_exact = None,
                   cmap = plt.cm.Greens,
                   label_ax = True,
                   cross_imm_labels = None,
                   **kwargs):
    fig = None
    
    if axis is None:
        fig, axis = plt.subplots()
    else:
        fig = axis.get_figure()

    if cross_imm_labels is None:
        cross_imm_labels = cross_immunities
        
    fmts = np.logspace(-10, 0, fineness)
    
    line_ids = np.linspace(crop_low, crop_high, len(cross_immunities))

    for c, line_id, label, in zip(cross_immunities,
                                  line_ids,
                                  cross_imm_labels):
        probs = [analytic_p_inocs(
            bottleneck,
            fmt,
            p_loss = 0,
            mut_neut_prob = wt_neut_prob * c,
            wt_neut_prob = wt_neut_prob,
            inoculum_size = inoculum_size,
            exact = inoculum_size < max_exact,
            verbose = False) /
                 analytic_p_inocs(bottleneck,
                                  fmt,
                                  p_loss = 0,
                                  inoculum_size=inoculum_size,
                                  wt_neut_prob=0,
                                  mut_neut_prob=0,
                                  exact = inoculum_size < max_exact)
                 for fmt in fmts]
        
        axis.plot(fmts, probs,
                  label=label,
                  lw=lw,
                  color=cmap(line_id),
                  **kwargs)

    if legend:
        leg = axis.legend()
        leg.set_title("cross-immunity")
        
    if label_ax:
        axis.set_xlabel("wild-type virion\n neutralization probability",
                        fontsize='xx-large')
        axis.set_ylabel("probability mutant present\nafter final bottleneck",
                        fontsize='xx-large')
        
    axis.set_title('mucus bottleneck $v$: {}\n'
                   'cell infection bottleneck $b$: {}'.format(
                       inoculum_size,
                       bottleneck),
                   loc='right', fontweight='light',
                   fontsize='medium')
    axis.set_xscale('log')
    axis.set_yscale('log')
