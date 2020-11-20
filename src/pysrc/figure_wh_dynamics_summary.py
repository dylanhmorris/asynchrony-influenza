#!/usr/bin/env python3

##################################################
# filename: figure_wh_dynamics_summary.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: generates multipanel summary
# of within-host dynamics results
##################################################

import sys

import plotting_style as ps
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
import parse_parameters as pp
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from matplotlib.ticker import ScalarFormatter

import pandas as pd

from figure_timing_heatmap import plot_heatmap
from plot_wh_timecourse import plot_wh_timecourse
from figure_wh_ngs import plot_wh_ngs, plot_ngs_hist
from plotting_functions import setup_multipanel, get_position, add_bounding_subplot



def plot_presentation_figures():
    non_trans_alpha = 0.4
    lineweight = 10
    mpl.rcParams['lines.linewidth'] = lineweight
    gen_param_file = "../../dat/RunParameters.mk"
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
    bottleneck=1

    fixed_results_dir = "../../out/within_host_results/minimal_fixed"
    var_results_dir = "../../out/within_host_results/minimal_visvar"
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
                  "E_w": False},
        "fixed": {"results-dir": fixed_results_dir,
                  "timecourse": "repl_visible",
                  "activation-time": 0,
                  "E_w": True},
        "var": {"results-dir": var_results_dir,
                "timecourse": "repl_visible",
                "activation-time": t_M,
                "E_w":True},
        "inoc": {"results-dir": var_results_dir,
                   "timecourse": "inoc_visible",
                   "activation-time": t_M,
                   "E_w": True}
    }
    

    for sim_type, metadata in sims_to_plot.items():
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        for ind, plot_freq in enumerate([False, True]):
            plot_axis = ax[ind]
            print("plotting {}...".format(sim_type))
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
                axis=plot_axis,
                t_M=metadata["activation-time"],
                E_w=metadata["E_w"],
                gen_param_file=gen_param_file,
                non_trans_alpha=non_trans_alpha,
                frequencies=plot_freq,
                analytical_frequencies=True)
            plot_axis.set_xlim([0, 8])
            plot_axis.set_xticks([0, 2, 4, 6, 8])
            if plot_freq:
                plot_axis.set_ylim([1e-6, 2])
        ax[0].set_ylabel("virions, cells")
        ax[1].set_ylabel("variant frequency")
        ax[0].set_xlabel("time (days)")
        ax[1].set_xlabel("time (days)")
        fig.tight_layout()
        plotname = 'timecourse_{}.pdf'.format(sim_type)
        fig.savefig('../../cons/fourth-year-talk/' + plotname)

    
def plot_presentation_legend():
    fig, ax = plt.subplots()
    ax.grid(b=False)
    non_trans_alpha = 0.4
    lineweight=5

    cells = mlines.Line2D([], [],
                          color=ps.cell_color,
                          lw=lineweight,
                          label='target\ncells')
    wt = mlines.Line2D([], [],
                          color=ps.wt_color,
                          lw=lineweight,
                          label='old variant\nvirus')
    mut = mlines.Line2D([], [],
                          color=ps.mut_color,
                          lw=lineweight,
                          label='new variant\nvirus')
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

    ax.legend(handles=handles,
              labels=labels,
              fontsize='x-large',
              frameon=True,
              ncol=1)


def main(output_path = None,
         fixed_results_dir = None,
         var_results_dir = None,
         gen_param_file = None,
         empirical_data_file = None,
         bottleneck = 1,
         heatmap_bottleneck = 1,
         increment = 0.05):
    
    if fixed_results_dir is None:
        fixed_results_dir = "../../out/within_host_results/minimal_visible"
    if var_results_dir is None:
        var_results_dir = "../../out/within_host_results/minimal_visvar"
    if output_path is None:
        output_path = "../../ms/main/figures/figure-wh-dynamics-summary.pdf"
    if gen_param_file is None:
        gen_param_file = "../../dat/RunParameters.mk"

    if empirical_data_file is None:
        empirical_data_file = "../../dat/cleaned/cleaned_wh_data.csv"

    aspect = .9
    width = 18.3
    height = width * aspect

    lineweight = width / 2.5
    mpl.rcParams['lines.linewidth'] = lineweight
    mpl.rcParams['font.size'] = width
    mpl.rcParams['legend.fontsize'] = "small"

    non_trans_alpha = 0.4

    params = pp.get_params(gen_param_file)
    detection_threshold = float(params["DEFAULT_DETECTION_THRESHOLD"])
    detection_limit = float(params["DEFAULT_DETECTION_LIMIT"])
    td_50 = float(params["DEFAULT_TD50"])
    transmission_cutoff = float(params["DEFAULT_TRANS_CUTOFF"])

    transmission_threshold = (
        -td_50 * np.log(1 - transmission_cutoff) /
        np.log(2))

    print(transmission_threshold)

    R0_wh = pp.get_param("R0_wh", "minimal_visvar", params)
    C_max = pp.get_param("C_max", "minimal_visvar", params)
    r_w = pp.get_param("r_w", "minimal_visvar", params)
    r_m = pp.get_param("r_m", "minimal_visvar", params)
    mu = pp.get_param("mu", "minimal_visvar", params)
    d_v = pp.get_param("d_v", "minimal_visvar", params)
    t_M = pp.get_param("t_M", "minimal_visvar", params)
    t_N = pp.get_param("t_N", "minimal_visvar", params)

    max_k = 8
    max_t_M = 3.5
    heatmap_t_final = 3

    # set up multi-panel figure
    fig = plt.figure(figsize=(width, height))

    n_cols = 4
    n_rows = 2
    height_ratios = [2, 3.5]
    ly2_inc = 0.07
    ly1 = 1 + 3 * ly2_inc / 2
    ly2 = 1 + ly2_inc
    lx = -0.05
    row_1_loc = (lx, ly1)
    row_2_loc = (lx, ly2)

    gs = gridspec.GridSpec(n_rows, n_cols,
                           height_ratios=height_ratios)
    all_heatmaps = add_bounding_subplot(fig,
                                        position=gs[0, 2:])

    all_plots = add_bounding_subplot(fig,
                                     position=gs[:,:])
    all_timecourses = add_bounding_subplot(fig,
                                           position=gs[1, :])

    
    plot_positions = [
        {"name": "empirical-detect",
        "grid_position": np.s_[0, 0],
        "sharex": None,
         "sharey": None,
         "letter_loc": row_1_loc},
        
        {"name": "empirical-hist",
         "grid_position": np.s_[0, 1],
         "sharex": None,
         "sharey": None,
         "letter_loc": row_1_loc},
        
        {"name": "one-percent-heatmap",
         "grid_position": np.s_[0, 2],
         "sharex": None,
         "sharey": None,
         "letter_loc": row_1_loc},
        
        {"name": "consensus-heatmap",
         "grid_position": np.s_[0, 3],
         "sharex": "one-percent-heatmap",
         "sharey": "one-percent-heatmap",
         "letter_loc": row_1_loc},
        
        {"name": "naive",
        "grid_position": np.s_[1, 0],
        "sharex": None,
         "sharey": None,
         "letter_loc": row_2_loc},
        
        {"name": "fixed",
         "grid_position": np.s_[1, 1],
         "sharex": None,
         "sharey": None,
         "letter_loc": row_2_loc},
        
        {"name": "var",
         "grid_position": np.s_[1, 2],
         "sharex": None,
         "sharey": None,
         "letter_loc": row_2_loc},
            
       {"name": "inoc",
        "grid_position": np.s_[1, 3],
        "sharex": None,
        "sharey": None,
         "letter_loc": row_2_loc}
                
    ]

    
    plots = setup_multipanel(fig,
                             plot_positions,
                             gridspec=gs)
    fig.tight_layout(rect=(0, 0.1, 1, 1),
                     w_pad=-0.5)

    for timecourse in ['naive', 'fixed', 'var', 'inoc']:
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1,
            hspace=0.15,
            subplot_spec=plots[timecourse])
        plots[timecourse + '-abs'] = plt.Subplot(fig, inner[0])
        plots[timecourse + '-freq'] = plt.Subplot(fig, inner[1],
                                                 sharex=plots.get('naive-freq', None),
                                                 sharey=plots.get('naive-freq', None))
        fig.add_subplot(plots[timecourse + '-abs'])
        fig.add_subplot(plots[timecourse + '-freq'])
        ax = plots[timecourse]
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.grid(b=False)
        ax.patch.set_alpha(0)
        ax.tick_params(labelcolor='w',
                       grid_alpha=0,
                       top=False,
                       bottom=False,
                       left=False,
                       right=False)
        ax.set_zorder(0)



    
    wt_col = "Vw"
    mut_col = "Vm"
    cell_col = "C"
    col_colors = [ps.wt_color, ps.mut_color, ps.cell_color]
    detection_linestyle = "dashed"
    detection_color="black"

    empirical_data = pd.read_csv(empirical_data_file)
    plot_wh_ngs(empirical_data,
                axis = plots['empirical-detect'],
                min_freq = 0.01,
                legend = True,
                edgecolor = "k")
    
    plot_ngs_hist(empirical_data,
                  axis = plots['empirical-hist'],
                  min_freq = 0,
                  legend = True)

    sims_to_plot = {
        "naive": {"results-dir": var_results_dir,
                  "timecourse": "naive",
                  "activation-time": t_N,
                  "E_w": False},
        "fixed": {"results-dir": fixed_results_dir,
                  "timecourse": "repl_visible",
                  "activation-time": 0,
                  "E_w": True},
        "var": {"results-dir": var_results_dir,
                "timecourse": "repl_visible",
                "activation-time": t_M,
                "E_w":True},
        "inoc": {"results-dir": var_results_dir,
                   "timecourse": "inoc_visible",
                   "activation-time": t_M,
                   "E_w": True}
    }
    
    for sim_type, metadata in sims_to_plot.items():
        for plot_freq, text in zip([False, True], ["abs", "freq"]):
            plot_axis_name = "{}-{}".format(
                sim_type, text)
            plot_axis = plots[plot_axis_name]
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
                axis=plot_axis,
                t_M=metadata["activation-time"],
                E_w=metadata["E_w"],
                gen_param_file=gen_param_file,
                non_trans_alpha=non_trans_alpha,
                frequencies=plot_freq,
                analytical_frequencies=True)
            plot_axis.set_xlim([0, 8])
            plot_axis.set_xticks([0, 2, 4, 6, 8])
            if plot_freq:
                plot_axis.set_ylim([1e-6, 2])

            ## remove inner labels
            if "naive" not in plot_axis_name:
                for label in plot_axis.get_yticklabels():
                    label.set_visible(False)
            if "abs" in plot_axis_name:
                for label in plot_axis.get_xticklabels():
                    label.set_visible(False)

    R0_wh = pp.get_param("R0_wh", "minimal_visvar", params) 
    d_v = pp.get_param("d_v", "minimal_visvar", params)

    min_prob = plot_heatmap(axis=plots["consensus-heatmap"],
                            contour_linewidth=(width/3) * 0.5,
                            contour_fontsize="large",
                            increment = increment,
                            mu = mu,
                            f_target = 0.5,
                            t_final = heatmap_t_final,
                            max_t_M = max_t_M,
                            max_k = max_k,
                            R0 = R0_wh,
                            c_w = 1,
                            c_m = 0,
                            d_v = d_v,
                            bottleneck=heatmap_bottleneck,
                            contour_levels=[1e-5, 1e-3, 1e-1],
                            cbar=True)
    
    plot_heatmap(axis=plots["one-percent-heatmap"],
                 contour_linewidth = (width/3) * 0.5,
                 contour_fontsize = "large",
                 increment = increment,
                 mu = mu,
                 bottleneck = heatmap_bottleneck,
                 f_target = 0.01,
                 min_prob = min_prob,
                 t_final = heatmap_t_final,
                 R0 = R0_wh,
                 d_v = d_v,
                 c_w = 1,
                 c_m = 0,
                 max_t_M = max_t_M,
                 max_k = max_k,
                 contour_levels = [1e-3, 1e-2, 1e-1],
                 cbar = False)

    star_size = 20
    star_marker = '*'
    star_facecolor = 'white'
    star_edgecolor = 'black'
    star_edgewidth = 1.5

    plots['empirical-detect'].set_ylabel('number of infections')
    plots['empirical-detect'].set_title('observed HA\npolymorphism')
    plots['empirical-hist'].set_title('variant within-host\nfrequencies')
    plots['empirical-hist'].set_xlabel('variant frequency')
    plots['empirical-hist'].set_ylabel('number of variants')

    for heatmap_name in ['consensus-heatmap',
                         'one-percent-heatmap']:
        hm = plots[heatmap_name]
        escape = 0.25
        star_x = 2.5
        sterilizing_k = (R0_wh - 1) * d_v
        fitness_diff = sterilizing_k * escape
        star_y = fitness_diff
        hm.plot(star_x, star_y,
                marker=star_marker,
                markersize=star_size,
                markerfacecolor=star_facecolor,
                markeredgewidth=star_edgewidth,
                markeredgecolor=star_edgecolor)
        hm.set_xlabel("")
        hm.grid(b=False)
        hm.set_yticks(np.arange(0, 8.5, 2))
        hm.set_xticks(np.arange(0, 3.5, 1))

        hm.set_xlabel("")

        ## need to reintroduce labels because
        ## axis sharing turns them off
        hm.xaxis.set_tick_params(
            labelbottom = True)

    ## need to reintroduce labels because
    ## axis sharing turns them off    
    plots['one-percent-heatmap'].yaxis.set_tick_params(
        labelleft = True)


    all_heatmaps.set_xlabel("time of recall response $t_M$ (days)")
    all_heatmaps.set_ylabel('selection strength $\delta$')

    
    all_timecourses.set_xlabel("time (days)", fontsize="xx-large")
    plots["naive-abs"].set_title("no recall\nresponse")
    plots["fixed-abs"].set_title("constant recall\nresponse")
    plots["var-abs"].set_title("recall response at\n48h")
    plots["inoc-abs"].set_title("response at 48h,\nvariant inoculated")
    plots["consensus-heatmap"].set_title("prob. new variant\n"
                                         "at consensus")
    plots["one-percent-heatmap"].set_title("prob. new variant\n"
                                           "at 1\%")
    
    plots["naive-abs"].set_ylabel("virions, cells")
    plots["naive-freq"].set_ylabel("variant frequency")

    plots['empirical-hist'].set_xlim([0, 0.5])
                
    # create legend
    cells = mlines.Line2D([], [],
                          color=ps.cell_color,
                          lw=lineweight,
                          label='target\ncells')
    wt = mlines.Line2D([], [],
                          color=ps.wt_color,
                          lw=lineweight,
                          label='old variant\nvirus')
    mut = mlines.Line2D([], [],
                          color=ps.mut_color,
                          lw=lineweight,
                          label='new variant\nvirus')
    ngs = mlines.Line2D([], [],
                          color='black',
                          lw=lineweight,
                          linestyle="dashed",
                          label='NGS detection\nlimit')
    analytical = mlines.Line2D([], [],
                               color='black',
                               lw=lineweight,
                               linestyle="dotted",
                               label='analytical new\nvariant frequency')
    antibody = mlines.Line2D([], [],
                             color=ps.immune_active_color,
                             lw=2*lineweight,
                             alpha=0.25,
                             label='recall response\nactive')
    star_leg = mlines.Line2D([], [],
                             color='white',
                             marker=star_marker,
                             markersize=star_size,
                             markerfacecolor=star_facecolor,
                             markeredgewidth=star_edgewidth,
                             markeredgecolor=star_edgecolor,
                             label="influenza-like\nparameters")
    
    handles = [cells, wt, mut, antibody, ngs, analytical, star_leg]
    labels = [h.get_label() for h in handles] 

    all_timecourses.legend(handles = handles,
                           labels = labels,
                           fontsize = 'x-large',
                           loc = "center",
                           bbox_to_anchor = (0.5, -.3),
                           frameon = False,
                           ncol = int(np.ceil(len(handles)/2)))    

    fig.savefig(output_path)



    
if __name__ == "__main__":
    if len(sys.argv) < 8:
        print("\nUSAGE: {} <fixed_results_dir> "
              "<var_results_dir> "
              "<empirical results file> "
              "<parameter file> "
              "<output path> <bottleneck> "
              "<heatmap n inoculated virions>"
              "\n"
              "".format(sys.argv[0]))
    else:
        main(fixed_results_dir = sys.argv[1],
             var_results_dir = sys.argv[2],
             empirical_data_file = sys.argv[3],
             gen_param_file = sys.argv[4],
             output_path = sys.argv[5],
             bottleneck = int(sys.argv[6]),
             heatmap_bottleneck = int(sys.argv[7]))
