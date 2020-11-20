#!/usr/bin/env python3

# makes heatmap of wh mutant frequency as a function of timing of adaptive response

from models.MinimalWithinHost import MinimalWithinHost
import models.MinimalWithinHost as mwh 
import numpy as np
import matplotlib as mpl
mpl.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from plotting_functions import plot_infection_course, plot_detectable_level, latex_pow_10, plot_df_timecourse
import plotting_style as ps
from plotting_style import wt_color, mut_color, cell_color, trans_threshold_color
import matplotlib.gridspec as gridspec
from distutils.sysconfig import parse_makefile
import pandas as pd
from wh_popgen import p_repl
import parse_parameters as pp

def texify(x, pos):
    return "{}".format(pos)




def get_repl_probs(
        ks,
        t_Ms,
        mu,
        f_target,
        t_final = None,
        c_w = None,
        c_m = None,
        bottleneck = None,
        R0 = None,
        d_v = None,
        verbose = False):
    n_ks = len(ks)
    n_t_Ms = len(t_Ms)
    
    results = np.zeros(n_ks * n_t_Ms).reshape((n_ks, -1))

    for i_k, k in enumerate(ks):
        if verbose:
            print(i_k)
        for j_t_M, t_M in enumerate(t_Ms):
            results[i_k, j_t_M] = p_repl(
                f_target,
                t_final,
                t_M,
                R0,
                d_v,
                bottleneck,
                k,
                c_w,
                c_m,
                mu)
    return results


def plot_heatmap(axis = None,
                 contour_linewidth = 5,
                 contour_fontsize = 20,
                 increment = 0.05,
                 mu = None,
                 f_target = None,
                 min_prob = 1,
                 t_final = None,
                 bottleneck = None,
                 R0 = None,
                 d_v = None,
                 c_w = None,
                 c_m = None,
                 clip_val=1e-10,
                 max_k = 8.5,
                 max_t_M = 3.5,
                 contour_levels = [1e-5, 1e-3, 1e-1],
                 cbar = True,
                 star = True,
                 verbose = False):
    if axis is None:
        fig, axis = plt.subplots()
    
    ks = np.arange(0, max_k, increment)
    t_Ms = np.arange(0, max_t_M, increment)
    origin_loc = "upper"

    a = get_repl_probs(ks,
                       t_Ms,
                       mu,
                       f_target,
                       bottleneck = bottleneck,
                       t_final = t_final,
                       R0 = R0,
                       d_v = d_v,
                       c_w = c_w,
                       c_m = c_m,
                       verbose = verbose)
    
    extent = [min(t_Ms) - increment/2, max(t_Ms) + increment/2,
              max(ks) + increment/2, min(ks) - increment/2]
    origin_loc = "upper"
    if a[a > clip_val].size > 0:
        a_min = np.min(a[a > clip_val])
    else:
        a_min = clip_val
    vmin_val = min(a_min, min_prob)
    heatmap = axis.imshow(a,
                          cmap='Reds',
                          norm = colors.LogNorm(
                              vmin = vmin_val,
                              vmax = 1),
                          alpha = 0.8,
                          extent = extent,
                          interpolation = 'nearest',
                          origin = origin_loc,
                          aspect = "auto")
    
    contours = axis.contour(a,
                            colors = "black",
                            levels = contour_levels,
                            linewidths = contour_linewidth,
                            extent = extent,
                            origin = origin_loc)

    clabs = axis.clabel(contours,
                        contours.levels,
                        inline = True,
                        fontsize = contour_fontsize,
                        fmt = LogFormatterMathtext(),
                        use_clabeltext = True)

    if cbar:
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right",
                                  size = "5%",
                                  pad = 0.05)
        cbar = plt.colorbar(heatmap,
                            cax = cax,
                            extend = 'both')
        cbar.ax.tick_params(labelsize = "large")

    axis.invert_yaxis()
    axis.set_xlabel("$t_M$")
    axis.set_ylabel("selection strength")
    axis.set_title("probability of selection\nto consensus")
    return vmin_val


def main(outdir,
         parameters = None):
    legend_fontsize = 30
    legend_handlelength = 0.3
    if parameters is not None:
        params = parse_makefile(parameters)

    # set plotting params
    axlabel_fontsize = 35
    ax_ticksize = 35
    t_final = 10
    
    width = 28
    height = 9

    # set up subplots
    gs = gridspec.GridSpec(1, 4, width_ratios=[1.08, 0.05, 1, 1])

    fig = plt.figure(figsize=(width, height))
    ax = [None, None, None]
    ax[0] = fig.add_subplot(gs[0])
    ax[1] = fig.add_subplot(gs[2])
    ax[2] = fig.add_subplot(gs[3], sharey=ax[1])

    # draw plots 
    plot_heatmap(axis=ax[0],
                 axlabel_fontsize=axlabel_fontsize,
                 contour_linewidth=5,
                 contour_fontsize=30)
    

    # style whole plot some more

    # no zeroth tick, for clash avoidance:

    y_max = max([ax[axis_no].get_ylim()[1] for axis_no in [1, 2]])
    for axis in [ax[1], ax[2]]:
        axis.set_xlim([0, t_final])
        axis.set_xlabel("Time (days)", fontsize=axlabel_fontsize)

    ax[1].set_ylabel("free virions")
    ax[1].set_ylabel("free virions")
    plt.setp(ax[2].get_yticklabels(), visible=False)
    for axis in ax:
        labels = axis.get_yticklabels()
        labels[0].set_visible(False)

    leg_labels = ["Uninfected cells",
                  "Wild-type virions",
                  "Mutant virions",
                  "Detection limit"]
    ax[2].legend(
        labels=leg_labels,
        fontsize=legend_fontsize,
        handlelength=legend_handlelength)

    fig.savefig(outdir)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("\nUSAGE: {} <output directory> <parameters>"
              "\n\n".format(sys.argv[0]))
    else:
        outdir = sys.argv[1]
        #parameters = sys.argv[2]
        main(outdir)
