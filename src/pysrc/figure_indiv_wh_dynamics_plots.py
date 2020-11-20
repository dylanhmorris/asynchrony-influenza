#!/usr/bin/env python3

##################################################
# filename: figure_indiv_wh_dynamics_plots.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: generates individual versions of
# the figures from figure_wh_dynamics_summary.py
# for use in conference talks etc
##################################################

import sys

import plotting_style as ps
import matplotlib as mpl
mpl.use("MacOSX")
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
import parse_parameters as pp
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines


from figure_timing_heatmap import plot_heatmap
from plot_wh_timecourse import plot_wh_timecourse
from plotting_functions import setup_multipanel, get_position

def main(output_path=None,
         fixed_results_dir=None,
         var_results_dir=None,
         gen_param_file=None,
         bottleneck=5,
         increment=0.5,
         landscape=True):
    if fixed_results_dir is None:
        fixed_results_dir = "../../out/within_host_results/minimal_visible"
    if var_results_dir is None:
        var_results_dir = "../../out/within_host_results/minimal_visvar"
    if output_path is None:
        output_path = "../../ms/main/figures/figure-wh-dynamics-summary.pdf"
    if gen_param_file is None:
        gen_param_file = "../../dat/RunParameters.mk"
    width = 18.3

    if landscape:
        height = 10
    else:
        height = 24.7

    heatmap_bottleneck = bottleneck
    lineweight = width / 2.5
    mpl.rcParams['lines.linewidth'] = lineweight
    non_trans_alpha = 0.4

    params = pp.get_params(gen_param_file)
    detection_threshold = float(params["DEFAULT_DETECTION_THRESHOLD"])
    transmission_threshold = float(params["DEFAULT_TRANS_THRESHOLD"])
    detection_limit = float(params["DEFAULT_DETECTION_LIMIT"])

    R0_wh = pp.get_param("R0_wh", "minimal_visvar", params)
    C_max = pp.get_param("C_max", "minimal_visvar", params)
    r_w = pp.get_param("r_w", "minimal_visvar", params)
    r_m = pp.get_param("r_m", "minimal_visvar", params)
    mu = pp.get_param("mu", "minimal_visvar", params)
    d_v = pp.get_param("d_v", "minimal_visvar", params)
    t_M = pp.get_param("t_M", "minimal_visvar", params)
    t_N = pp.get_param("t_N", "minimal_visvar", params)
    
    n_pairs = 5

    wt_col = "Vw"
    mut_col = "Vm"
    cell_col = "C"
    col_colors = [ps.wt_color, ps.mut_color, ps.cell_color]
    detection_linestyle = "dashed"
    detection_color="black"

    sims_to_plot = {
        "naive": {"results-dir": var_results_dir,
                  "timecourse": "naive",
                  "activation-time": t_N,
                  "title": "no recall\nresponse",
                  "E_w": False},
        "fixed": {"results-dir": fixed_results_dir,
                  "timecourse": "repl_visible",
                  "title": "constant recall\nresponse",
                  "activation-time": 0,
                  "E_w": True},
        "var": {"results-dir": var_results_dir,
                "timecourse": "repl_visible",
                "title": "delayed recall\nresponse",
                "activation-time": t_M,
                "E_w":True},
        "inoc": {"results-dir": var_results_dir,
                   "timecourse": "inoc_visible",
                   "activation-time": t_M,
                   "title": "delayed, mutant\ninoculated",
                   "E_w": True}
    }
    
    for sim_type, metadata in sims_to_plot.items():
        fig, plot_axes = plt.subplots(1, 2, figsize=(16, 8))
        for plot_freq, text in zip([0, 1], ["abs", "freq"]):
            plot_wh_timecourse(
                metadata["timecourse"],
                bottleneck,
                metadata["results-dir"],
                wt_col,
                mut_col,
                col_colors=col_colors,
                cell_col=cell_col,
                detection_threshold=detection_threshold,
                detection_color=detection_color,
                detection_limit=detection_limit,
                detection_linestyle=detection_linestyle,
                transmission_threshold=transmission_threshold,
                axis=plot_axes[plot_freq],
                t_M=metadata["activation-time"],
                E_w=metadata["E_w"],
                gen_param_file=gen_param_file,
                non_trans_alpha=non_trans_alpha,
                frequencies=plot_freq,
                analytical_frequencies=True)
            plot_axes[plot_freq].set_xlim([0, 8])
            plot_axes[plot_freq].set_xticks([0, 2, 4, 6, 8])
        plot_axes[0].set_ylim(ymax=1e10)
        plot_axes[1].set_ylim(ymax=5)

        plot_axes[0].set_ylabel("virions, cells")
        plot_axes[0].set_xlabel("time (days)")
        plot_axes[1].set_xlabel("time (days)")
        plot_axes[1].set_ylabel("variant frequency")
        fig.tight_layout()
        fig.savefig("../../out/wh-plot-{}.pdf"
                    "".format(sim_type))
    heatmap_fig, heatmaps = plt.subplots(1, 2, figsize=(16, 8))
    min_prob = plot_heatmap(axis=heatmaps[0],
                            axlabel_fontsize=20,
                            contour_linewidth=(width/3) * 0.75,
                            contour_fontsize="x-large",
                            increment=increment,
                            mu=mu,
                            f_target=0.5,
                            bottleneck=heatmap_bottleneck,
                            params=params,
                            cbar=True)
    
    plot_heatmap(axis=heatmaps[1],
                 axlabel_fontsize=20,
                 contour_linewidth=(width/3) * 0.75,
                 contour_fontsize="x-large",
                 increment=increment,
                 mu=mu,
                 bottleneck=heatmap_bottleneck,
                 f_target=0.01,
                 min_prob=min_prob,
                 params=params,
                 cbar=True)
    heatmaps[0].set_xlabel("time of recall response")
    heatmaps[1].set_xlabel("time of recall response")

    heatmaps[0].set_ylabel("selection strength")
    heatmaps[1].set_title("probability of selection\n to 1%")
    heatmap_fig.tight_layout()
    heatmap_fig.savefig("../../out/talk-heatmaps.pdf")

def make_legend():
    width = 18.3
    lineweight = width / 2.5

    leg_fig = plt.figure()
    cells = mlines.Line2D([], [],
                          color=ps.cell_color,
                          lw=lineweight,
                          label='target\ncells')
    wt = mlines.Line2D([], [],
                          color=ps.wt_color,
                          lw=lineweight,
                          label='wild-type\nvirus')
    mut = mlines.Line2D([], [],
                          color=ps.mut_color,
                          lw=lineweight,
                          label='mutant\nvirus')
    ngs = mlines.Line2D([], [],
                          color='black',
                          lw=lineweight,
                          linestyle="dashed",
                          label='NGS detection\nlimit')
    analytical = mlines.Line2D([], [],
                               color='black',
                               lw=lineweight,
                               linestyle="dotted",
                               label='analytical mutant\nfrequency')
    antibody = mlines.Line2D([], [],
                             color=ps.immune_active_color,
                             lw=2*lineweight,
                             alpha=0.25,
                             label='recall response\nactive')
    handles = [cells, wt, mut, antibody, ngs, analytical]
    labels = [h.get_label() for h in handles] 
    plt.legend(handles=handles,
               labels=labels,
               loc="center",
               frameon=True)
    leg_fig.savefig("../../out/talk-wh-legend.pdf")
