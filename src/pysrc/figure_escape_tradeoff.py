#!/usr/bin/env python3

##############################################
# name: figure_escape_tradeoff.py
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

from point_of_transmission import analytic_p_inocs, neutralize_prob_from_z, poisson_max_sus

from plotting_functions import add_bounding_subplot

mpl.rcParams['font.size'] = 25

def plot_escape_tradeoff(axis = None,
                         bottleneck = None,
                         f_mut = None,
                         inoculum_size = None,
                         cross_immunities = [None],
                         legend = True,
                         lw = 5,
                         plot_drift = True,
                         crop_low = 0.4,
                         crop_high = 0.8,
                         fineness = 50,
                         max_exact = 1,
                         cmap = plt.cm.Greens,
                         label_ax = True,
                         cross_imm_labels = None,
                         inoculum_model = None,
                         **kwargs):
    fig = None
    
    if axis is None:
        fig, axis = plt.subplots()
    else:
        fig = axis.get_figure()

    if cross_imm_labels is None:
        cross_imm_labels = cross_immunities

    wt_neut_probs = np.linspace(0, 1, fineness)
            
    line_ids = np.linspace(crop_low, crop_high, len(cross_immunities))

    for c, line_id, label, in zip(cross_immunities,
                                  line_ids,
                                  cross_imm_labels):
        mut_neut_probs = [
            c * wt_neut_prob for wt_neut_prob in wt_neut_probs]

        probs = [analytic_p_inocs(
            bottleneck,
            f_mut,
            p_loss = 0,
            mut_neut_prob = mut_neut_prob,
            wt_neut_prob = wt_neut_prob,
            inoculum_size = inoculum_size,
            exact = inoculum_size < max_exact,
            verbose = False,
            inoculum_model = inoculum_model)
                 for wt_neut_prob, mut_neut_prob in zip(
                 wt_neut_probs, mut_neut_probs)]
        
        axis.plot(wt_neut_probs,
                  probs,
                  label = label,
                  lw = lw,
                  color = cmap(line_id),
                  **kwargs)

    if plot_drift:
        axis.axhline(analytic_p_inocs(bottleneck,
                                      f_mut,
                                      p_loss = 0,
                                      inoculum_size = inoculum_size,
                                      wt_neut_prob = 0,
                                      mut_neut_prob = 0,
                                      inoculum_model = inoculum_model),
                     color = ps.drift_color,
                     linestyle = 'dotted',
                     lw = lw/2)

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

    
def main(inoculum_mults = [1, 10, 25],
         bottleneck_sizes = [1, 5, 10],
         outpath = '../../ms/main/figures/figure-escape-tradeoff.pdf',
         param_path = None,
         fineness = 200,
         implication_z_m = None,
         inoculum_model = "poisson"):
    nrows = len(bottleneck_sizes)
    ncols = len(inoculum_mults)
    width = 22
    height = (nrows / ncols) * width
    
    if param_path is None:
        param_path = "../../dat/RunParameters.mk"
    params = pp.get_params(param_path)

    f_mut = pp.get_param("f_mut", "", params)
    print("fmt: ", f_mut)

    fig, axes = plt.subplots(nrows=nrows,
                             ncols=ncols,
                             sharex=True,
                             sharey= True, # 'row',
                             figsize=(width, height))

    full_plot = add_bounding_subplot(fig)
    xticklabs = ["" for tick in full_plot.get_xticks()]
    yticklabs = ["" for tick in full_plot.get_yticks()]
    full_plot.set_xticklabels(xticklabs)
    full_plot.set_yticklabels(yticklabs)
    full_plot.tick_params(pad=70)
    
    prob_ticks = [0, 0.25, 0.5, 0.75, 1]
    
    for i, bn in enumerate(bottleneck_sizes):
        for j, mult in enumerate(inoculum_mults):
            print("plotting final bottleneck b: {}, v: {}..."
                  "".format(bn, bn * mult))
            v = bn * mult

            plot_escape_tradeoff(
                f_mut = f_mut,
                axis = axes[i, j],
                bottleneck = bn,
                inoculum_size = bn * mult,
                cross_immunities = [0, 0.5, 0.9, 0.99],
                legend = False,
                label_ax = False,
                cross_imm_labels = [
                    "0\%",
                    "50\%",
                    "90\%",
                    "99\%"],
                fineness = fineness,
                inoculum_model = inoculum_model)
            axes[i, j].label_outer()
            axes[i, j].set_xticks(prob_ticks)
            axes[i, j].set_yscale('log')
            axes[i, j].set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
    leg = axes[0, 0].legend(ncol = 2,
                            handlelength = 1)
    leg.set_title("sIgA cross immunity\n($\sigma = \kappa_m / \kappa_w$)")
    full_plot.set_ylabel("probability mutant present\n"
                         "after final bottleneck")
    full_plot.set_xlabel("wild-type neutralization "
                         "probability $\kappa_{w}$")
    fig.tight_layout()
    fig.savefig(outpath)

    
if __name__ == "__main__":
    usage_msg = ("USAGE: {} <parameter file> <output path>"
                 "\n\n".format(sys.argv[0]))
    if len(sys.argv) < 3:
        print(usage_msg)
    elif not os.path.isfile(sys.argv[1]):
        print(usage_msg)
    else:
        main(param_path=sys.argv[1],
             outpath=sys.argv[2])
