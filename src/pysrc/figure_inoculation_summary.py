#!/usr/bin/env python3

##############################################
# name: figure_inoculation_summary.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: summary figure containing
# various comparisons of the probability
# of replication selection and inoculation
# selection as a function of various parameters
##############################################

import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines

import numpy as np
from scipy.stats import norm

import plotting_style as ps
import analyze_selection_comparison as asc
import parse_parameters as pp
import point_of_transmission as pot

from plotting_functions import setup_multipanel, get_position, add_bounding_subplot, subdivided_hist, plot_image_from_path

from figure_mucosal_filter import plot_survive_bottleneck
from figure_distribution_shift import plot_filtered_inocula, plot_sim_inocula

from figure_optimal_selector import plot_optimal_selector_cluster, plot_sus_funcs
from figure_cutdown import plot_cutdown


def setup_figure():
    magnif = 3
    width = magnif * 18.3
    nrows = 3
    ncols = 4
    height = width
    font_mag = 1
    
    lineweight = width / 3
    mpl.rcParams['lines.linewidth'] = lineweight
    mpl.rcParams['axes.formatter.limits'] = (-3, 6)
    mpl.rcParams['font.size'] *= magnif * font_mag
    mpl.rcParams['ytick.major.pad'] =  width / 4
    mpl.rcParams['xtick.major.pad'] = width / 2

    # set up multi-panel figure
    fig = plt.figure(figsize=(width, height))
    gs = gridspec.GridSpec(nrows, ncols)

    
    schematic_row, inoc_dist_row, selectors_row = 0, 1, 2
    
    prefilter_col, postfilter_col, mut_inf_col = 0, 1, 2

    virion_survival_col, cutdown_col, sus_mod_col, opt_sel_col = 0, 1, 2, 3

    plot_positions = [
        {"name": "bottleneck-schematic",
         "grid_position": np.s_[schematic_row, :],
         "letter_loc": (-0.0125, 1.1),
         "sharex": None,
         "sharey": None},

        {"name": "virion-survival",
         "grid_position": np.s_[selectors_row, virion_survival_col],
         "sharex": None,
         "sharey": None},
        
        {"name": "cutdown",
         "grid_position": np.s_[selectors_row, cutdown_col],
         "sharex": None,
         "sharey": None},

        {"name": "susceptibility-models",
         "grid_position": np.s_[selectors_row, sus_mod_col],
         "sharex": None,
         "sharey": None},

        {"name": "optimal-selectors",
         "grid_position": np.s_[selectors_row, opt_sel_col],
         "sharex": "susceptibility-models",
         "sharey": None}
    ]
    
    plots = setup_multipanel(fig,
                             plot_positions,
                             gridspec = gs,
                             verbose = True,
                             letters = ['a', 'e', 'f', 'g', 'h', 'i'])

    distrow = gs[inoc_dist_row, :]
    inner_gs = gridspec.GridSpecFromSubplotSpec(
        1, 3,
        subplot_spec=distrow)
    plots['all_dists'] = add_bounding_subplot(fig,
                                              position=distrow)

    letters = ['B', 'C', 'D']
    for idist, dist in enumerate(
            ['pre-filter', 'post-filter', 'mutant-creator']):
        plotname = dist + '-dist'
        plots[plotname] = fig.add_subplot(
            inner_gs[idist],
            sharex=plots.get('pre-filter-dist', None),
            sharey=plots.get('pre-filter-dist', None))
        letter_x, letter_y = ps.letter_loc
        plots[plotname].text(letter_x, letter_y, letters[idist],
                             transform = plots[plotname].transAxes,
                             fontsize = ps.letter_size,
                             fontweight = 'bold',
                             va='top')
        

    return (fig, plots)


def plot_inoculation_summary(
        wh_sim_data_dir=None,
        schematic_path=None,
        run_parameter_path=None,
        file_output_path=None,
        escape = 0.75,
        v = 10,
        b = 1,
        sigma = 0.75,
        inoculum_model = 'poisson',
        z_homotypic = None):
    
    width = 18.3
    lineweight = width / 3
    span_fontsize = 'x-large'
    schem_fontsize = 'xx-large'

    fig, plots = setup_figure()

    # get within-host data
    dataframe = asc.get_vb_data(wh_sim_data_dir)

    # get model name
    model_name = os.path.basename(wh_sim_data_dir)
    print("model name: ", model_name)
    
    # get parameters
    params = pp.get_params(run_parameter_path)
    
    inoculum_sizes = [1, 3, 10, 50, 100, 200]

    ## plot bottleneck schematic
    plot_image_from_path(schematic_path,
                         axis = plots['bottleneck-schematic'],
                         retain_axis = True)
    

    schem_labs = [
        "",
        "excretion\nbottleneck",
        "inter-host\nbottleneck",
        "mucus\nbottleneck",
        "sIgA\nbottleneck",
        "cell infection\nbottleneck"]
    
    plots['bottleneck-schematic'].set_xticklabels(
        schem_labs,
        fontsize = schem_fontsize)

    schem_ticks = plots['bottleneck-schematic'].get_xticks()

    schem_ticks = [(x + schem_ticks[i + 1]) / 2
                   if x < len(schem_ticks) - 1 else
                   x + (x - schem_ticks[i - 1]) / 2
                   for i, x in enumerate(schem_ticks)]
    schem_ticks = schem_ticks[:len(schem_labs)]
    print(schem_ticks)
    plots['bottleneck-schematic'].set_xticks(
        schem_ticks)
    plots['bottleneck-schematic'].set_xlim(
        left = (schem_ticks[0] + schem_ticks[1]) / 2)
    plots['bottleneck-schematic'].grid(b = False)

    
    plots['bottleneck-schematic'].set_yticklabels([])

    ## get parameters and plot filtered inocula
    z_wt = pp.get_param("z_wt",
                        model_name,
                        params)

    if z_homotypic is None:
        z_homotypic = z_wt

    f_mut = pp.get_param("f_mut", "", params)
    print("fmt: ", f_mut)

    mut_wt_neut_ratio = pp.get_param("mut_wt_neut_ratio",
                                     model_name,
                                     params)
    
    kappa_ws = pot.neutralize_prob_from_z(z_wt,
                                          np.array(inoculum_sizes),
                                          "poisson")
    sim_p_mut = dataframe.groupby('n_encounters_iga')['p_mut_inoc'].mean()[inoculum_sizes]
    plot_filtered_inocula(
        inoculum_sizes,
        sim_p_mut,
        axis=plots["pre-filter-dist"])
    plots['pre-filter-dist'].legend(
        ['neither', 'old variant', 'new variant', 'both'],
        loc='lower right')

    plot_filtered_inocula(
        inoculum_sizes,
        sim_p_mut,
        axis=plots["post-filter-dist"],
        kappa_w = kappa_ws,
        kappa_m = kappa_ws * mut_wt_neut_ratio)
    print("sim_p_mut:", sim_p_mut)

    plot_sim_inocula(
        dataframe,
        inoculum_sizes=inoculum_sizes,
        axis=plots["mutant-creator-dist"],
        f_mut_threshold=0.5)


    print("plotting virion survival...")
                    
    plot_survive_bottleneck(
        axis = plots["virion-survival"],
        f_mut = f_mut,
        bottleneck = b,
        lw = lineweight * 2.5,
        drift = True,
        inoculum_scaleup = v / b,
        cmap = plt.cm.Greens,
        fineness = 10)
    plots['virion-survival'].set_ylim(bottom = 0)
    plots['virion-survival'].set_xlim(left = 0,
                                      right = 1)

    
    print("Plotting susceptibility model examples...")
    sus_func_cmaps = [
        plt.cm.Purples,
        plt.cm.Greens]

    sus_line_alpha = 0.8

    plot_sus_funcs(
        escape = escape,
        cmaps = sus_func_cmaps,
        axis = plots["susceptibility-models"],
        z_homotypic = z_homotypic,
        line_alpha = sus_line_alpha)
    
    print("Plotting optimal selector by memory age...")
    c_darks = [0.7, 0.4]
    drift_style = "dashed"
    
    plot_optimal_selector_cluster(
        f_mut = f_mut,
        inoculum_size = v,
        bottleneck = b,
        escape = escape,
        axis = plots["optimal-selectors"],
        cmaps = sus_func_cmaps,
        darkness = c_darks[0],
        line_alpha = sus_line_alpha,
        inoculum_model = inoculum_model,
        z_homotypic = z_homotypic,
        drift_style = drift_style,
        legend = False)

    plot_optimal_selector_cluster(
        f_mut = f_mut,
        inoculum_size = v,
        bottleneck = b,
        escape = escape,
        axis = plots["optimal-selectors"],
        cross_imm_sigma = sigma,
        cmaps = sus_func_cmaps,
        darkness = c_darks[1],
        line_alpha = sus_line_alpha,
        inoculum_model = inoculum_model,
        z_homotypic = z_homotypic,
        plot_drift = False,
        legend = False)

    # manual legend for optimal selector panel
    sel_cluster_handles = [
        mlines.Line2D([], [],
                      color = plt.cm.Greys(c_darks[0]),
                      label = "implied\nby $z_{m}$"),
        mlines.Line2D([], [],
                      color = plt.cm.Greys(c_darks[1]),
                      label = ("{} ".format(sigma) +
                               "$\kappa_{w}$"))
    ]
    plots['optimal-selectors'].legend(handles = sel_cluster_handles,
                                      title = "$\kappa_m$",
                                      fancybox = True,
                                      frameon = True)
    
    print("plotting cutdown at iga bottleneck...")
    plot_cutdown(axis=plots['cutdown'], f_mut = f_mut)
    plots['cutdown'].set_xticks([1, 100, 200])


    
    # plot styling
    cluster_ticks = [0, 1, 2, 3, 4, 5]
    for plotname in ["susceptibility-models",
                     "optimal-selectors"]:
        plots[plotname].set_xticks(cluster_ticks)
        plots[plotname].set_xlim(left = cluster_ticks[0],
                                 right = cluster_ticks[-1])
        plots[plotname].set_ylim(bottom = 0)
        plots[plotname].set_xlabel(
            "distance between\n"
            "host memory and old variant")

    for plotname in ["pre-filter-dist",
                     "post-filter-dist",
                     "mutant-creator-dist"]:
        plots[plotname].label_outer()
        plots[plotname].set_yticks([0, 0.25, 0.5, 0.75, 1])

    plots['pre-filter-dist'].set_ylabel('proportion of inocula')
    plots['pre-filter-dist'].set_title('before IgA\n')
    plots['post-filter-dist'].set_title('after IgA\n')
    plots['mutant-creator-dist'].set_title('new variant emerged\n')

    emergence_y_text = "new variant infections\nper inoculation"
    plots["optimal-selectors"].set_ylabel(
        emergence_y_text)
    
    plots["virion-survival"].set_ylabel(
        "prob. new variant survives bottleneck")

    plots["virion-survival"].set_ylabel(
        "prob. new variant survives bottleneck")

    plots["all_dists"].set_xlabel('virions encountering IgA\n($v$)',
                                  fontsize = span_fontsize)

    

    fig.tight_layout(h_pad=-0.1)
    print("saving figure to {}".format(file_output_path))
    fig.savefig(file_output_path)

if __name__=="__main__":
    if len(sys.argv) < 5:
        print("USAGE: {} <schematic path> "
              "<wh sim data dir> "
              "<param file> <output path> "
              "\n\n".format(sys.argv[0]))
    else:
        plot_inoculation_summary(
            schematic_path=sys.argv[1],
            wh_sim_data_dir=sys.argv[2],
            run_parameter_path=sys.argv[3],
            file_output_path=sys.argv[4])
