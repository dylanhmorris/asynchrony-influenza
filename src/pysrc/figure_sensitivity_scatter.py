#!/usr/bin/env python3

###################################################
# filename: figure_sensitivty_scatter.py
# author: Dylan H. Morris <dhmorris@princeton.edu>
# 
# description: generates extended data figure
# showing sensitivity analysis scatterplot in
# varied parameters
###################################################

import sys
import os

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import ScalarFormatter

import analyze_selection_comparison as asc
import  plotting_functions as pf
import plotting_style as ps
import parse_parameters as pp

mpl.rcParams['lines.markeredgewidth'] = 1

order_of_outcomes = ['mut-inf',
                     'ratio']

palette = [ps.repl_color,
           ps.inocs_color]

observation_alpha = 0.8

## denominator for frequencies

denom = 10000

def calc_quantities(x):
    n_detectable = np.sum(x.transmissible)
    n_new_variant = np.sum((x.transmissible) &
                           (x.peak_mutant_freq >= 0.01))
    n_repl = np.sum((x.transmissible) &
                    (x.peak_mutant_freq >= 0.01) &
                    ~(x.is_inoc))
    
    n_inoc = np.sum((x.transmissible) &
                    (x.peak_mutant_freq >= 0.01) &
                    (x.is_inoc))
    n_total = x.size

    
    return (
        pd.Series(
            [
                n_new_variant * 100 / n_detectable,

                (n_inoc >= n_repl) * 1
            ],
        
            index=order_of_outcomes)
        )

def get_tidy_data(results_dir):
    dat = asc.get_bottleneck_data(results_dir)
    parms = asc.get_sensitivity_paramsets(results_dir)
    parms['paramset_id'] = parms['paramset_id'].astype('int')
    ## extract needed simulation info and get to tidy format for plotting

    individual_results = dat.groupby(
        ['bottleneck',
         'paramset_id']).apply(calc_quantities)

    grouped_results = dat.groupby(
        'bottleneck').apply(calc_quantities)

    individual_results['niter'] = dat.groupby(
        ['bottleneck',
         'paramset_id'])['rng_seed'].count()
    
    individual_results = pd.melt(individual_results.reset_index(),
                                 id_vars=['bottleneck',
                                          'paramset_id',
                                          'niter',
                                          'ratio'],
                                 var_name='outcome',
                                 value_name='frequency')
    return(
        pd.merge(individual_results, parms,
                 on=['paramset_id', 'bottleneck']))

def parm_scatter(tidy_results,
                 parm,
                 outcome,
                 axis=None,
                 cmap=palette,
                 observation_alpha=observation_alpha,
                 size=50,
                 xlims=None):

    if axis is None:
        fig, axis = plt.subplots()
        
    df_ind = tidy_results[
        (tidy_results.outcome == outcome)]

    axis.scatter(
        df_ind[parm],
        df_ind['frequency'],
        alpha=observation_alpha,
        color=[cmap[int(x)] for x in tidy_results['ratio']],
        linewidth=0.8,
        s=size,
        marker='o',
        edgecolors='k')

    if xlims is None:
        xlims = [df_ind[parm].min(), df_ind[parm].max()]

    min_val, max_val = xlims[0], xlims[1]
    x_range = max_val - min_val
    xlim_vals = [min_val - 0.1 * x_range,
                 max_val + 0.1 * x_range]        
    axis.set_xlim(xlim_vals)


def parm_striplot(tidy_results,
                  parm,
                  outcome,
                  axis=None,
                  cmap=palette,
                  observation_alpha=observation_alpha,
                  size=7,
                  jitter=True,
                  xlims=None,
                  legend=False):

    if axis is None:
        fig, axis = plt.subplots()
        
    df_ind = tidy_results[
        (tidy_results.outcome == outcome)]

    # hack to make sure that absent observations
    # don't lead to something getting the wrong hue
    pruned_cmap = [cmap[int(x)] for x in np.sort(pd.unique(df_ind['ratio']))]
    sorted_order = np.sort(pd.unique(
        tidy_results[parm]))
    
    sns.stripplot(
        x=parm,
        y="frequency",
        hue="ratio",
        ax=axis,
        data=df_ind,
        order=sorted_order,
        alpha=observation_alpha,
        palette=pruned_cmap,
        jitter=jitter,
        dodge=True,
        linewidth=0.8,
        size=size,
        edgecolor='k')

    if not legend:
        axis.get_legend().remove()


def main(results_dir=None,
         param_path=None,
         output_path=None):

    width = 18.3
    aspect = 1
    height = width * aspect
    figsize = (width, height)
    legend_markersize = width
    
    if results_dir is None:
        results_dir = "../../out/sensitivity_analysis_results/minimal_visvar"

    if param_path is None:
        param_path = "../../dat/RunParameters.mk"
        
    if output_path is None:
        output_path = ("../../ms/supp/"
                       "figs/figure-sensitivity-scatter-visvar.pdf")

    parm_dict = pp.get_params(param_path)
    model_name = os.path.basename(results_dir)
    
    indiv = get_tidy_data(results_dir)
    extra_columns = ['paramset_id', 'niter', 'outcome', 'frequency', 'ratio']
    parms = [col for col in indiv.columns if col not in extra_columns]
    n_parms = len(parms)
    
    n_cols = 4
    n_rows = int(np.ceil(n_parms / n_cols))
    n_axes = n_rows * n_cols

    fig, ax = plt.subplots(nrows=n_rows, ncols=n_cols,
                           figsize=figsize,
                           sharey=True)
    for ind, parm in enumerate(parms):
        axis = ax.flatten()[ind]

        if parm != 'bottleneck':
            xlim = [
                pp.get_param(param_name = parm + "_MIN",
                             model_name = model_name,
                             param_dict = parm_dict),
                pp.get_param(param_name = parm + "_MAX",
                             model_name = model_name,
                             param_dict = parm_dict)]
            parm_scatter(indiv, parm, 'mut-inf', axis=axis,
                         xlims=xlim)
            
        else:
            xlim = None
            parm_striplot(indiv, parm, 'mut-inf',
                          axis=axis, xlims=xlim,
                          jitter=0.2)
            
        axis.set_xlabel(ps.parm_display_names.get(parm,
                                                  parm))
        axis.set_yscale('symlog')
        axis.set_ylim([-0.1, 110])

    # legend in an unused axis
    leg_ax = ax.flatten()[n_parms]
    labels = ['de novo new variants\nmore common',
              'inoculated new variants\nmore common']

    legend_elements = [mpl.lines.Line2D(
        [0], [0],
        marker='o',
        color=color,
        alpha=observation_alpha,
        lw=0,
        markersize=legend_markersize,
        markeredgecolor='k',
        label=lab) for color, lab in
                       zip(palette, labels)]
    leg_ax.legend(handles=legend_elements)
    leg_ax.spines['top'].set_color('none')
    leg_ax.spines['bottom'].set_color('none')
    leg_ax.spines['left'].set_color('none')
    leg_ax.spines['right'].set_color('none')
    leg_ax.grid(b=False)
    leg_ax.tick_params(
        labelcolor='w',
        grid_alpha=0,
        top=False,
        bottom=False,
        left=False,
        right=False)

    # delete unused axes
    for i_ax in range(n_parms + 1, n_axes):
        fig.delaxes(ax.flatten()[i_ax])

    # style axes
    for axis in ax[:, 0].flatten():
        axis.set_ylabel('new variant infections\n'
                        'per 100 detectable infections')
        axis.yaxis.set_major_formatter(ScalarFormatter())
    for axis in ax[:, 1:].flatten():
        plt.setp(axis.get_yticklabels(), visible=False)
        
    fig.tight_layout()
    fig.savefig(output_path)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("\nUSAGE: {} <sensitivity results> "
              "<run params> "
              "<output path>"
              "\n"
              "".format(sys.argv[0]))
    else:
        main(sys.argv[1],
             sys.argv[2],
             sys.argv[3])
