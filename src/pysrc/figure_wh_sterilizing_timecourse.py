#!/usr/bin/env python3

###################################################
# filename: figure_wh_sterilizing_timecourse
# author: Dylan H. Morris <dhmorris@princeton.edu>
# 
# description: generates extended data figure
# showing what sterilizing neutralization during
# the timecourse of infection looks like 
###################################################

import sys

import plotting_style as ps
import plotting_functions as pf
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import parse_parameters as pp
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import pandas as pd

import analyze_selection_comparison as asc
import point_of_transmission as pot
import wh_popgen as whp

from plot_wh_timecourse import plot_wh_timecourse

analytic_linestyle = 'dashed'
prob_alpha = 0.5

denom_exponent = 5 # number of infections per for symlog
denom = 10 ** denom_exponent

aspect = .5
width = 18.3
height = width * aspect

lineweight = width / 2.5
mpl.rcParams['lines.linewidth'] = lineweight
mpl.rcParams['font.size'] = 0.8 * width
mpl.rcParams['axes.titlesize'] = "medium"
mpl.rcParams['legend.fontsize'] = "medium"


possible_outcomes = [
    'no detectable\ninfection',
    'detectable infection with\nde novo new variant',
    'detectable infection with\ninoculated new variant',
    'detectable infection with\nwild-type']


def calc_apply(x):
    denom = np.sum(x.transmissible | (~x.transmissible))
    return pd.Series(
        [np.sum(~x.transmissible) / denom,
         np.sum((x.transmissible) &
                 (x.peak_mutant_freq >= 0.01) &
                 (~x.is_inoc)) / denom,
         np.sum((x.transmissible) &
                 (x.peak_mutant_freq >= 0.01) &
                 (x.is_inoc)) / denom,
         np.sum((x.transmissible) &
                 (x.peak_mutant_freq < 0.01)) / denom
        ],
        index=possible_outcomes)


def calc_mut_inf(x):
    n_detectable = np.sum(x.transmissible)
    n_new_variant = np.sum((x.transmissible) &
                           (x.peak_mutant_freq >= 0.01))
    return n_new_variant / n_detectable

def calc_stoch_loss(x):
    n_mut = np.sum(x.ever_mutated & (~x.is_inoc))
    n_surv = np.sum((x.emergence_time < 99) & (~x.is_inoc))
    return n_surv / n_mut

def calc_mutated(x):
    return np.sum(x.ever_mutated & (~x.is_inoc))

def calc_raw_repl(x):
    return np.sum(x.transmissible & (~x.is_inoc))

def calc_repl(x):
    return np.sum(x.transmissible & (~x.is_inoc)) / np.sum(
        x.transmissible | (~x.transmissible))

def calc_p_inocs(f_mut, bn, R0_wh):
    p_inoculated = (1 - (1 - f_mut) ** bn)
    p_loss = pot.minimal_p_loss(f_mut, R0_wh, bn)
    return p_inoculated * (1 - p_loss)


def setup_figure():

    # set up multi-panel figure
    fig = plt.figure(figsize=(width, height))

    n_cols = 5
    n_rows = 2
    letter_y = 1.15
    letter_x = -0.05

    gs = gridspec.GridSpec(n_rows, n_cols)
    
    all_plots = pf.add_bounding_subplot(
        fig,
        position=gs[:,:])

    all_timecourses = pf.add_bounding_subplot(
        fig,
        position=gs[:, :3])

    
    plot_positions = [
        
        {"name": "sterilized-abs",
         "grid_position": np.s_[0, 0],
         "sharex": None,
         "sharey": None},
        
        {"name": "replication-selected-abs",
         "grid_position": np.s_[0, 1],
         "sharex": "sterilized-abs",
         "sharey": "sterilized-abs"},
        
        {"name": "inoculation-selected-abs",
         "grid_position": np.s_[0, 2],
         "sharex": "sterilized-abs",
         "sharey": "sterilized-abs"},

        {"name": "proportion-plot",
         "grid_position": np.s_[:, 3:],
         "sharex": None,
         "sharey": None,
         "letter_loc": (letter_x, letter_y + (1 - letter_y) / 2),
         "include_letter": True},
        
        {"name": "sterilized-freq",
         "grid_position": np.s_[1, 0],
         "sharex": "sterilized-abs",
         "sharey": None,
         "include_letter": False},
        
        {"name": "replication-selected-freq",
         "grid_position": np.s_[1, 1],
         "sharex": "sterilized-abs",
         "sharey": "sterilized-freq",
         "include_letter": False},
        
        {"name": "inoculation-selected-freq",
         "grid_position": np.s_[1, 2], 
         "sharex": "sterilized-abs",
         "sharey": "sterilized-freq",
         "include_letter": False}
    ]

    
    plots = pf.setup_multipanel(fig,
                                plot_positions,
                                gridspec = gs,
                                letter_loc = (letter_x, letter_y))

    plots['all-plots'] = all_plots
    plots['all-timecourses'] = all_timecourses
    
    fig.tight_layout()
    return (fig, plots)


def main(output_path = None,
         results_dir = None,
         gen_param_file = None,
         bottleneck = 200):
    if results_dir is None:
        results_dir = "../../out/within_host_results/minimal_sterilizing"
    if output_path is None:
        output_path = "../../ms/main/figures/figure-wh-sterilizing-timecourse.pdf"
    if gen_param_file is None:
        gen_param_file = "../../dat/RunParameters.mk"


    fig, plots = setup_figure()

    ## read in parameters
    params = pp.get_params(gen_param_file)
    detection_threshold = float(params["DEFAULT_DETECTION_THRESHOLD"])
    transmission_threshold = float(params["DEFAULT_TRANS_THRESHOLD"])
    detection_limit = float(params["DEFAULT_DETECTION_LIMIT"])
    f_mut_default = float(params["DEFAULT_F_MUT"])

    model_name = 'minimal_sterilizing'
    
    R0_wh = pp.get_param("R0_wh", model_name, params)
    C_max = pp.get_param("C_max", model_name, params)
    r_w = pp.get_param("r_w", model_name, params)
    r_m = pp.get_param("r_m", model_name, params)
    mu = pp.get_param("mu", model_name, params)
    d_v = pp.get_param("d_v", model_name, params)
    k = pp.get_param("k", model_name, params)
    cross_imm = pp.get_param("cross_imm", model_name, params)
    t_M = pp.get_param("t_M", model_name, params)
    t_N = pp.get_param("t_N", model_name, params)
    p_loss = pot.minimal_p_loss(3e-5, R0_wh, 1)
    print("p loss:", p_loss)

    # set styling

    wt_col = "Vw"
    mut_col = "Vm"
    cell_col = "C"
    col_colors = [ps.wt_color, ps.mut_color, ps.cell_color]
    detection_linestyle = "dashed"
    detection_color="black"
    non_trans_alpha = 0.4

    sims_to_plot = {'sterilized': 'repl_invisible',
                    'replication-selected': 'repl_visible',
                    'inoculation-selected': 'inoc_visible'}
    
    for sim_type, timecourse in sims_to_plot.items():
        for plot_freq, text in zip([False, True], ["abs", "freq"]):
            plot_axis_name = "{}-{}".format(
                sim_type, text)
            plot_axis = plots[plot_axis_name]

            print("plotting {}...".format(sim_type))
            plot_wh_timecourse(
                timecourse,
                bottleneck,
                results_dir,
                wt_col,
                mut_col,
                col_colors = col_colors,
                cell_col = cell_col,
                detection_threshold = detection_threshold,
                detection_color = detection_color,
                detection_limit = detection_limit,
                detection_linestyle = detection_linestyle,
                transmission_threshold = transmission_threshold,
                axis = plot_axis,
                t_M = t_M,
                E_w = True,
                gen_param_file = gen_param_file,
                non_trans_alpha = non_trans_alpha,
                frequencies = plot_freq,
                analytical_frequencies = True)
            plot_axis.set_xlim([0, 8])
            plot_axis.set_xticks([0, 2, 4, 6, 8])
            if plot_freq:
                plot_axis.set_ylim([1e-6, 2])
            if ("sterilized" not in plot_axis_name):
                for label in plots[plot_axis_name].get_yticklabels():
                    label.set_visible(False)
            
            if "abs" in plot_axis_name:
                for label in plots[plot_axis_name].get_xticklabels():
                    label.set_visible(False)


    dat = asc.get_bottleneck_data(results_dir)
    f_muts = asc.get_mean_fmut(dat)
    bottlenecks = np.sort(pd.unique(dat.bottleneck))

    analytic_p_inocs = np.array([
        calc_p_inocs(float(f_muts.at[bn]),
                     bn,
                     R0_wh)
        for bn in bottlenecks
    ])
        
    analytic_p_repls = np.array([
        whp.p_repl_declining(
            bn,
            mu,
            R0_wh,
            d_v,
            k,
            1,
            cross_imm)
        for bn in bottlenecks
    ])


    analytic_p_repls = analytic_p_repls * (1 - analytic_p_inocs)
    
    analytic_p_wts = np.zeros_like(analytic_p_repls)
    
    analytic_p_elses = 1 - analytic_p_inocs - analytic_p_repls

    ## extract needed simulation info
    ## and get to tidy format for plotting

    print('analytical p repl:')
    print(analytic_p_repls)

    print('raw simulated repl events:')
    print(dat.groupby('bottleneck').apply(calc_raw_repl))
    
    print('simulated repl probs:')
    print(dat.groupby('bottleneck').apply(calc_repl))

    results = dat.groupby('bottleneck').apply(calc_apply)
    
    results = pd.melt(results.reset_index(), id_vars='bottleneck',
                      var_name='outcome',
                      value_name='frequency')

    my_palette = ['grey', ps.repl_color, ps.inocs_color,
                  ps.wt_color]

    analyticals = [
        analytic_p_elses,
        analytic_p_repls,
        analytic_p_inocs,
        analytic_p_wts
    ]

    for outcome, color, analytical in zip(possible_outcomes,
                                          my_palette,
                                          analyticals):
        df = results[results.outcome == outcome]
        plots['proportion-plot'].plot(
            df.bottleneck,
            df.frequency * denom,
            marker = 'o',
            markersize = 20,
            linestyle = "",
            markeredgecolor = 'k',
            color = color,
            alpha = prob_alpha,
            label = outcome)
        
        plots['proportion-plot'].plot(
            bottlenecks,
            analytical * denom,
            marker = '+',
            markeredgewidth = 5,
            markersize = 20,
            linestyle = "",
            alpha = 1,
            color = color)


    # style resulting plots
    plots['proportion-plot'].set_yscale('symlog')
    plots['proportion-plot'].set_title('distribution of outcomes\n')
    plots['proportion-plot'].legend(loc='center left',
                                    bbox_to_anchor=(0, 0.7))
    plots['proportion-plot'].set_xscale('log')
    plots['proportion-plot'].set_xlabel('bottleneck')
    plots['proportion-plot'].set_ylabel('frequency per $10^{}$ inoculations'
                                        ''.format(denom_exponent))

    plots['sterilized-abs'].set_ylabel('virions, cells')
    plots['sterilized-freq'].set_ylabel('frequency')

    plots['sterilized-abs'].set_title(
        'no detectable infection\n')

    plots['replication-selected-abs'].set_title(
        'detectable infection\nwith de novo new variant')
    
    plots['inoculation-selected-abs'].set_title(
        'detectable infection\nwith inoculated new variant')
    
    plots['all-timecourses'].set_xlabel('time (days)')
    plots['proportion-plot'].set_xticks(bottlenecks)
    plots['proportion-plot'].set_xticklabels(bottlenecks)

    fig.tight_layout()
    fig.savefig(output_path)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("\nUSAGE: {} <results dir> <parameter file>"
              "<output path>\n"
              "".format(sys.argv[0]))
    else:
        main(results_dir=sys.argv[1],
             gen_param_file=sys.argv[2],
             output_path=sys.argv[3])
