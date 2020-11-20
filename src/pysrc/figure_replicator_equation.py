#!/usr/bin/env python3

######################################################
# name: figure_replicator_equation.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: plot trajectories for
# the replicator equation
######################################################

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import wh_popgen as whp
import plotting_style as ps
import plotting_functions as pf


def main(argv):
    if len(argv) < 2:
        print("USAGE: ./{} <outpath>\n\n"
              "outpath: path to output figure\n\n")
        return 0

    outpath = argv[1]
    
    width = 12
    height = width
    color_low, color_high = 0.3, 0.7
    lineweight = width / 2.5
    mpl.rcParams['lines.linewidth'] = lineweight
    mpl.rcParams['font.size'] = 1.5 * width
    mpl.rcParams['legend.fontsize'] = "x-small"
    mpl.rcParams['legend.loc'] = "upper left"

    fig = plt.figure(figsize = (width, height))

    plot_positions = [
        {"name": "varying-selection",
         "row": 1,
         "col": 1,
         "sharex": None,
         "sharey": None},
        
        {"name": "varying-init-freq",
         "row": 1,
         "col": 2,
         "sharex": "varying-selection",
         "sharey": "varying-selection"},
        
        {"name": "varying-selection-ongoing",
         "row": 2,
         "col": 1,
         "sharex": "varying-selection",
         "sharey": None},

        {"name": "varying-init-freq-ongoing",
         "row": 2,
         "col": 2,
         "sharex": "varying-selection-ongoing",
         "sharey": "varying-selection-ongoing"},

    ]

    all_plots = pf.add_bounding_subplot(fig,
                                        position = "111")
    plots = pf.setup_multipanel(fig,
                                plot_positions)

    times = np.linspace(0, 4, 1000)

    deltas = [12, 9, 6, 3]
    delta_colors = np.linspace(color_high,
                               color_low,
                               len(deltas))
    delta_cmap = plt.cm.Greens

    mu = 0.33e-5
    R0 = 5
    d_v = 4
    t_peak = 1.5
    t_M = 2
    t_emerge = 0
    
    mu_mult_exps = [3, 2, 1, 0]
    mu_colors = np.linspace(color_high,
                            color_low,
                            len(mu_mult_exps))
    mu_cmap = plt.cm.Blues


    for delta, color in zip(deltas, delta_colors):
        freqs = [whp.freq_mut_select(t, mu, delta)
                 for t in times]
        freqs_ongoing = [whp.freq_mut_ongoing_piecewise(
            t,
            t_emerge,
            t_M,
            t_peak,
            mu,
            delta,
            mu,
            R0 * d_v) for t in times]
        
        plots["varying-selection"].plot(
            times,
            freqs,
            color = delta_cmap(color),
            label = "${}$".format(delta))
        plots["varying-selection-ongoing"].plot(
            times,
            freqs_ongoing,
            color = delta_cmap(color))
        plots["varying-selection-ongoing"].axvline(
            t_M,
            color = "k",
            linestyle = "dashed")


    delta_fix = 6
    
    for mult_exp, color in zip(mu_mult_exps, mu_colors):
        freqs = [whp.freq_mut_select(t,
                                     mu * 10**(mult_exp),
                                     delta_fix)
                 for t in times]
        freqs_ongoing = [whp.freq_mut_ongoing_piecewise(
            t,
            t_emerge,
            t_M,
            t_peak,
            mu * 10**(mult_exp),
            delta_fix,
            mu,
            R0 * d_v) for t in times]

        plots["varying-init-freq"].plot(
            times,
            freqs,
            color = mu_cmap(color),
            label = "$\\mu \\times 10^{}$".format(mult_exp))
        plots["varying-init-freq-ongoing"].plot(
            times,
            freqs_ongoing,
            color = mu_cmap(color))
        plots["varying-init-freq-ongoing"].axvline(
            t_M,
            color = "k",
            linestyle = "dashed")

        
    plots["varying-selection"].legend(
        title = "$\delta$",
        handlelength = width / 10)
    plots["varying-init-freq"].legend(
        title = "$f_0$",
        handlelength = width / 20)

    freq_ticks = [0, 0.25, 0.5, 0.75, 1]
    for plotname, plot in plots.items():
        plot.label_outer()
        plot.set_xlim([0, 4])
        if not "ongoing" in plotname:
            plot.set_yticks(freq_ticks)
            plot.set_yticklabels(["${}$".format(x)
                                  for x in freq_ticks])
        else:
            plot.set_yscale("log")
            
    all_plots.set_xlabel("time since first new variant (days)")
    plots["varying-selection"].set_ylabel("variant frequency $f_m$")
    plots["varying-selection-ongoing"].set_ylabel("variant frequency $f_m$")
    
    fig.tight_layout()

    fig.savefig(outpath)

    return 0


if __name__ == "__main__":
    main(sys.argv)
