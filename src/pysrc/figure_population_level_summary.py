#!/usr/bin/env python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import matplotlib.lines as mlines
import numpy as np

import plotting_style as ps
import parse_parameters as pp
from plotting_functions import setup_multipanel

from point_of_transmission import neutralize_prob_from_z, epidemic_p_inocs, inocs_by_pct_immune
from multi_class_final_size import multi_class_inoculations
from models.susceptibility_models import pick_sus_func

def plot_epidemic_p_inocs(axis = None,
                          cmap = plt.cm.Reds,
                          distribution = [],
                          f_mut = None,
                          bottleneck = None,
                          v = None,
                          population_R0s = None,
                          escape = None,
                          susceptibility_model = "linear",
                          N = 1,
                          cross_imm_sigma = None,
                          z_wt = None,
                          legend = True,
                          leg_cmap = None,
                          **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    wt_pcts = np.arange(0, 1.0, 0.01)
    
    line_ids = np.linspace(0.15, 0.7, len(population_R0s))

    # setup for legend
    if leg_cmap is None:
        leg_cmap = cmap
        
    handlelines = []

    for R0_population, line_id in zip(population_R0s, line_ids):
        p_inocs_vals = [epidemic_p_inocs(
            wt_pct,
            f_mut = f_mut,
            bottleneck = bottleneck,
            R0_population = R0_population,
            escape = escape,
            distribution = distribution,
            susceptibility_model = susceptibility_model,
            v = v,
            z_wt = z_wt,
            cross_imm_sigma = cross_imm_sigma,
            N = N)
                        for wt_pct in wt_pcts]
        
        axis.plot(wt_pcts, p_inocs_vals,
                  color=cmap(line_id),
                  **kwargs)
        handlelines.append(
            mlines.Line2D([], [],
                          color = leg_cmap(line_id),
                          label = str(R0_population),
                          **kwargs))

    axis.set_xlabel("level of immunity")
    axis.set_ylabel("probability of\nnew variant chain")
    axis.set_xlim([0, 1])
    axis.set_ylim(bottom = 0)
    if legend:
        legend = axis.legend(handles = handlelines,
                             frameon = True,
                             fancybox = True,
                             framealpha = 1,
                             labelspacing = 0.25,
                             handlelength = 0.5)
        legend.set_title("$\mathfrak{R}_0$")
    axis.grid(b = True)


def plot_inoculation_selection_vs_immunity(
        axis = None,
        f_mut = None,
        wt_neut_prob = None,
        mut_neut_probs = [None],
        final_bottleneck = None,
        mucus_bottleneck = None,
        inoculum_model = None,
        bottleneck_model = None,
        legend = True,
        x_min = 0,
        x_max = 1,
        fineness = 1000,
        cmap = plt.cm.Reds,
        leg_cmap = None,
        col_low = 0.4,
        col_high = 0.8,
        **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    if leg_cmap is None:
        leg_cmap = cmap

    n_mut_neut_probs = len(mut_neut_probs)
    colors = np.linspace(col_low, col_high, n_mut_neut_probs)
    for line_id, mut_neut_prob in enumerate(mut_neut_probs):
        xs, ys = inocs_by_pct_immune(
            f_mut = f_mut,
            final_bottleneck = final_bottleneck,
            mucus_bottleneck = mucus_bottleneck,
            wt_neut_prob = wt_neut_prob,
            mut_neut_prob = mut_neut_prob,
            x_min = x_min,
            x_max = x_max,
            fineness = fineness,
            inoculum_model = inoculum_model,
            bottleneck_model = bottleneck_model)

        axis.plot(xs, ys, label = mut_neut_prob,
                  color=cmap(colors[line_id]),
                  **kwargs)

    if legend:
        leg_lines = []
        for line_id, mut_neut_prob in enumerate(mut_neut_probs):
            leg_lines += [Line2D([0], [0],
                                 color=leg_cmap(colors[line_id]),
                                 **kwargs)]

        axis.legend(leg_lines, mut_neut_probs,
                    handlelength=0.5,
                    title="$\kappa_m$")


def plot_reinoculations(axis=None,
                        cmap=plt.cm.Reds,
                        leg_cmap=None,
                        legend=True,
                        R0s = [1.2, 1.5, 1.8],
                        **kwargs):
    if not axis:
        fig, axis = plt.subplots()

    handlelines = []

    if leg_cmap is None:
        leg_cmap = cmap

    s0s = np.arange(0, 1.01, 0.01)[::-1]
    line_ids = np.linspace(0.15, 0.7, len(R0s))

    lineweight = ps.standard_lineweight
    axis_fontsize = ps.standard_axis_fontsize
    axis_ticksize = ps.standard_axis_ticksize

    for R0, line_id in zip(R0s, line_ids):
        xs = 1 - s0s
        reinoculations = [multi_class_inoculations([1, 0],
                                                   [1 - x, x],
                                                   R0)[1][1]
                          for x in xs]
        axis.plot(xs, reinoculations,
                  color=cmap(line_id),
                  **kwargs)

        handlelines.append(
            mlines.Line2D([], [],
                          color=leg_cmap(line_id),
                          label=str(R0),
                          **kwargs))

    axis.set_xlabel("level of immunity")
    axis.set_ylabel("reinoculations per capita")
    axis.set_xlim([0, 1])
    axis.set_ylim(bottom=0)

    if legend:
        legend = axis.legend(handles=handlelines,
                             frameon=True,
                             fancybox=True,
                             framealpha=1,
                             labelspacing=0.25,
                             handlelength=0.5)
        legend.set_title("$\mathfrak{R}_0$")
        
    axis.grid(b=True)

    
def main(param_file=None,
         output_path=None):

    ## figure styling / setup
    width = 18.3
    mpl.rcParams['font.size'] = width / 1.5
    mpl.rcParams['lines.linewidth'] = width / 2.5
    kernel_lw = width / 4.5
    mpl.rcParams['ytick.major.pad']= width / 2
    mpl.rcParams['xtick.major.pad']= width / 2
    mpl.rcParams['legend.fontsize']= width * 0.9
    mpl.rcParams['legend.title_fontsize']= width 
    mpl.rcParams['legend.handlelength']= 2

    height = width / 3
    fig = plt.figure(figsize=(width, height))

    nrows = 1
    ncols = 3
    gs = gridspec.GridSpec(nrows, ncols)

    ## multipanel setup
    pop_row = 0
    e_rate_col, reinoc_col, epi_escape_col = 0, 1, 2
    
    plot_positions = [
        {"name": "emergence-rate",
         "grid_position": np.s_[pop_row, e_rate_col],
         "sharex": None,
         "sharey": None},
        
        {"name": "reinoculations",
         "grid_position": np.s_[pop_row, reinoc_col],
         "sharex": None,
         "sharey": None},
        
        {"name": "epi-escape",
         "grid_position": np.s_[pop_row, epi_escape_col],
         "sharex": None,
         "sharey": None}
    ]

    letter_loc = (-0.1, 1.15)
    plots = setup_multipanel(fig,
                             plot_positions,
                             letter_loc=letter_loc,
                             gridspec=gs)



    ## parametrization
    params = pp.get_params(param_file)
    b = int(params.get("DEFAULT_BOTTLENECK"))
    v = int(float(params.get("DEFAULT_VB_RATIO")) * b)
    f_mut = float(params.get("DEFAULT_F_MUT"))
    z_wt = float(params.get("DEFAULT_Z_WT"))
    k = float(params.get("CONSTANT_CHAIN_K"))
    mu = float(params.get("DEFAULT_MU"))
    d_v = float(params.get("DEFAULT_D_V"))
    R0 = float(params.get("DEFAULT_R0_WH"))

    v_mult = 5
    epi_sigma = 0.75
    escape = 1
    
    print("parsed parameter file: v = {}, b = {}".format(v, b))
    print("parsed parameter file: f_mut = {}".format(f_mut))
    print("parsed parameter file: z_wt = {}".format(z_wt))
    
    pop_R0s = [1.5, 2, 2.5]
    print("Plotting reinoculations vs population immunity...")
    plot_reinoculations(
        axis = plots["reinoculations"],
        R0s = pop_R0s,
        cmap = plt.cm.Greens,
        leg_cmap = plt.cm.Greys)

    print("Plotting new chains vs population immunity...")
    plot_epidemic_p_inocs(
        axis = plots["epi-escape"],
        f_mut = f_mut,
        bottleneck = b,
        cross_imm_sigma = epi_sigma,
        v = v,
        population_R0s = pop_R0s,
        escape = escape,
        leg_cmap = plt.cm.Greys,
        z_wt = z_wt,
        legend = True)

    print("Plotting new chains vs population immunity with v = {}..."
          "".format(v * v_mult))

    plot_epidemic_p_inocs(
        axis = plots["epi-escape"],
        cross_imm_sigma = epi_sigma,
        f_mut = f_mut,
        bottleneck = b,
        v = v * v_mult,
        population_R0s = pop_R0s,
        cmap = plt.cm.Purples,
        linestyle = 'dotted',
        escape = escape,
        z_wt = z_wt,
        legend = False)    
    
    print("Plotting emergence vs population immunity...")
    mut_neut_prob_list = [0.75, 0.9, 0.99]
    
    plot_inoculation_selection_vs_immunity(
        axis = plots["emergence-rate"],
        f_mut = f_mut,
        mut_neut_probs = mut_neut_prob_list,
        wt_neut_prob = neutralize_prob_from_z(
            z_wt,
            v,
            "poisson"),
        final_bottleneck = b,
        mucus_bottleneck = v,
        cmap = plt.cm.Reds,
        linestyle = "solid",
        legend = True,
        leg_cmap=plt.cm.Greys,
        inoculum_model = "poisson",
        bottleneck_model = "binomial")

    plot_inoculation_selection_vs_immunity(
        axis = plots["emergence-rate"],
        f_mut = f_mut,
        mut_neut_probs = mut_neut_prob_list,
        wt_neut_prob = neutralize_prob_from_z(
            z_wt,
            v * v_mult,
            "poisson"),
        final_bottleneck = b,
        mucus_bottleneck = v * v_mult,
        cmap = plt.cm.Purples,
        linestyle = "dotted",
        legend = False,
        inoculum_model = "poisson",
        bottleneck_model = "binomial")

    #####################################
    # plot styling
    #####################################
    imm_level_text = "population fraction immune\nto wild-type"
    emergence_y_text = "new variant infections\nper inoculation"
    frac_ticks = [0, 0.25, 0.5, 0.75, 1]
    frac_tick_labs = ["$0$", "$0.25$", "$0.5$", "$0.75$", "$1$"]
    
    plots["reinoculations"].set_xlabel(
        imm_level_text)
    plots["reinoculations"].set_ylabel(
        "per capita reinoculations")
    plots["epi-escape"].set_xlabel(
        imm_level_text)
    plots["epi-escape"].set_ylabel(
        "per capita probability\n"
        "of new variant infection")
    plots["emergence-rate"].set_xlabel(
        imm_level_text)
    plots["emergence-rate"].set_ylabel(
        emergence_y_text)
    plots["epi-escape"].set_ylim([0, 3e-4])


    for plot in plots.values():
        plot.set_xlim(left = 0, right = 1)
        plot.set_ylim(bottom = 0)
        plot.set_xticks(frac_ticks)
        plot.set_xticklabels(frac_tick_labs)


    
        
    for plot in plots.values():
        if not type(plot.yaxis._scale) == mpl.scale.LogScale:
            plot.set_ylim(bottom=0)
    
    plt.tight_layout()

    # save
    fig.savefig(output_path)

    return 0


if __name__ == "__main__":
    if len(sys.argv) < 1 + 2:
        print("USAGE: {} "
              "<param file> <output path> "
              "\n\n".format(sys.argv[0]))
    else:
        main(param_file=sys.argv[1],
             output_path=sys.argv[2])


