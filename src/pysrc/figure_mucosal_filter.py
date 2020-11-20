
#!/usr/bin/env python3

##############################################
# name: figure_mucosal_filter.py
# author: Dylan Morris <dhmorris@princeton.edu>
#
##############################################

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import binom
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plotting_style as ps
import parse_parameters as pp
from plotting_functions import setup_multipanel

from point_of_transmission import analytic_p_inocs

def plot_survive_bottleneck(axis=None,
                            crop_low=0.4,
                            crop_high=0.8,
                            cmap=plt.cm.Greens,
                            wt_neut_probs=[0.25, 0.5, 0.75, 1],
                            f_mut=3e-5,
                            bottleneck=1,
                            inoculum_scaleup=1,
                            min_imm=0,
                            drift=True,
                            linestyles = ['solid'] * 4,
                            drift_cmap=plt.cm.Greys,
                            fineness=50,
                            lw=2,
                            **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    line_ids = np.linspace(crop_low, crop_high, len(wt_neut_probs))
    mut_neut_prob = np.linspace(0, 1, fineness)
    
    print("plotting bottleneck {}".format(bottleneck))
    for line_id, line_style, wt_neut_prob in zip(line_ids,
                                                 linestyles,
                                                 wt_neut_probs):
        inoculum_size=inoculum_scaleup * bottleneck
        y_vals = [analytic_p_inocs(bottleneck,
                                   f_mut,
                                   p_loss = 0,
                                   inoculum_size=inoculum_size,
                                   wt_neut_prob=wt_neut_prob,
                                   mut_neut_prob=mut_neut_prob)
                  for mut_neut_prob in mut_neut_prob]
        axis.plot(mut_neut_prob,
                  y_vals,
                  color=cmap(line_id),
                  label=wt_neut_prob,
                  lw=lw,
                  linestyle=line_style,
                  **kwargs)

        if drift:
            axis.axhline(analytic_p_inocs(bottleneck,
                                          f_mut,
                                          p_loss = 0,
                                          inoculum_size=inoculum_size,
                                          wt_neut_prob=0,
                                          mut_neut_prob=0),
                         color=drift_cmap(line_id),
                         linestyle='dashed',
                         lw=lw)
        
    leg = axis.legend(frameon = True,
                      fancybox = True,
                      handlelength = lw / 10)
    leg.set_title("$\kappa_w$")
    axis.set_ylabel("prob. mutant survives bottleneck")
    axis.set_xlabel("mutant neutralization\nprob. ($\kappa_m$)")

def presentation_figure(outdir="../../out/mucosal-filter.pdf"):
    width = 18.3
    height = 6
    lw = width / 2.5

    fig, ax= plt.subplots(1, 2, figsize=(width, height))

        
    plot_survive_bottleneck(axis=ax[0],
                            lw=lw)
    plot_survive_bottleneck(axis=ax[1],
                            lw=lw,
                            inoculum_scaleup=10)
    ax[0].set_title("$v=1$, $b=1$")
    ax[1].set_title("$v=10$, $b=1$")
    plt.tight_layout()
    fig.savefig(outdir)
