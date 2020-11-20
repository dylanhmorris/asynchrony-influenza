#!/usr/bin/env python3

###################################################
# filename: figure_sensitivty.py
# author: Dylan H. Morris <dhmorris@princeton.edu>
# 
# description: generates extended data figure
# showing sensitivity analysis results on p_inocs
# and p_repl
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

mpl.rcParams['lines.markeredgewidth'] = 1

order_of_outcomes = ['detectable infection',
                     'new variant infection',
                     'new variant per detectable',
                     'ratio']

palette = [ps.repl_color,
           ps.inocs_color]

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
                np.mean((x.transmissible)),

                n_new_variant / n_total,
                
                n_new_variant / n_detectable,
                
                (n_inoc >= n_repl) * 1
            ],
        
            index=order_of_outcomes)
        )

def get_tidy_data(results_dir):
    dat = asc.get_bottleneck_data(results_dir)
    parms = asc.get_sensitivity_paramsets(results_dir)
    parms['paramset_id'] = parms['paramset_id'].astype('int')
    ## extract needed simulation info and get to tidy format for plotting
    dat = pd.merge(dat, parms, on=['paramset_id', 'bottleneck'])

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
    
    grouped_results = pd.melt(grouped_results.reset_index(),
                              id_vars=['bottleneck'],
                              var_name='outcome',
                              value_name='frequency')

    return (individual_results, grouped_results, parms)


def results_plot(tidy_results,
                 outcome,
                 palette=None,
                 order=order_of_outcomes,
                 axis=None,
                 observation_alpha=0.5,
                 jitter=0.2,
                 legend=False,
                 legend_bbox_to_anchor=None):

    if axis is None:
        fig, axis = plt.subplots()

    # if axes is None:
    #     fig, axes = plt.subplots(ncols=outcomes.size)
    # elif len(axes) != outcomes.size:
    #     raise ValueError("Number of axes must be same "
    #                      "as number of unique results "
    #                      "in tidy data")

    sorted_bottlenecks = np.sort(
        pd.unique(tidy_results.bottleneck))
        
    df_ind = tidy_results[
        tidy_results.outcome == outcome]

    # hack to make sure that absent observations
    # don't lead to something getting the wrong hue
    pruned_palette = palette
        
    sns.stripplot(
        x="bottleneck",
        y="frequency",
        hue="ratio",
        ax=axis,
        data=df_ind,
        order=sorted_bottlenecks,
        alpha=observation_alpha,
        palette=pruned_palette,
        jitter=jitter,
        dodge=True,
        linewidth=1.2,
        size=15,
        edgecolor='k')
    
    #axis.set_yscale('log')
    handles, labels = axis.get_legend_handles_labels()

    if legend:
        axis.legend(handles[0:],
                    ['more de novo\nnew variants',
                     'more inoculated\nnew variants'],
                    markerscale=3,
                    bbox_to_anchor=legend_bbox_to_anchor)
    else:
        axis.get_legend().remove()


def setup_figure(n_rows = 2,
                 n_cols = 3):
    aspect = 0.8
    width = 18.3
    height = width * aspect

    lineweight = width / 2.5
    mpl.rcParams['lines.linewidth'] = lineweight
    mpl.rcParams['font.size'] = width * 0.9
    mpl.rcParams['axes.titlesize'] = 'large'
    mpl.rcParams['axes.titleweight'] = 'bold'

    # set up multi-panel figure
    fig = plt.figure(figsize=(width, height))

    gs = gridspec.GridSpec(n_rows, n_cols)
    
    all_plots = pf.add_bounding_subplot(
        fig,
        position=gs[:,:])

    sensitivity_model_names = ['fixed', 'visvar']

    row_plots = [
        pf.add_bounding_subplot(
            fig,
            position = gs[row, :])
        for row in range(n_rows)
    ]
    
    plot_positions = [
        
        {"name": "fixed-any-inf",
         "grid_position": np.s_[0, 0],
         "sharex": None,
         "sharey": None},
        
        {"name": "visvar-any-inf",
         "grid_position": np.s_[1, 0],
         "sharex": "fixed-any-inf",
         "sharey": "fixed-any-inf"},
            

        {"name": "fixed-mut-inf",
         "grid_position": np.s_[0, 1],
         "sharex": "fixed-any-inf",
         "sharey": None},
        
        {"name": "visvar-mut-inf",
         "grid_position": np.s_[1, 1],
         "sharex": "fixed-any-inf",
         "sharey": "fixed-mut-inf"},

        {"name": "fixed-mut-inf-per",
         "grid_position": np.s_[0, 2],
         "sharex": "fixed-any-inf",
         "sharey": None},
        
        {"name": "visvar-mut-inf-per",
         "grid_position": np.s_[1, 2],
         "sharex": "fixed-any-inf",
         "sharey": "fixed-mut-inf-per"}
    ]
    
    plots = pf.setup_multipanel(fig,
                                plot_positions,
                                gridspec=gs)

    for row, name in enumerate(sensitivity_model_names):
        plots[name] = row_plots[row]

    plots['all'] = all_plots
    
    fig.tight_layout()
    return (fig, plots)


def main(results_prefix,
         output_path):

    if results_prefix is None:
        results_prefix = "../../out/sensitivity_analysis_results"

    ## setup plot 
    fig, plots = setup_figure()

    ## do plotting 
    for ind, paramset_type in enumerate(['fixed', 'visvar']):
        results_dir = os.path.join(results_prefix, 'minimal_' + paramset_type)
        indiv, overall, parms = get_tidy_data(results_dir)

        suffixes = [
            '-any-inf',
            '-mut-inf',
            '-mut-inf-per'
        ]

        outcomes_to_plot = [
            'detectable infection',
            'new variant infection',
            'new variant per detectable']


        palettes = [
            [ps.wt_color],
            
            [ps.repl_color,
             ps.inocs_color],
            
            [ps.repl_color,
             ps.inocs_color],
        ]
             
        for outcome, suffix, palette in zip(outcomes_to_plot,
                                            suffixes,
                                            palettes):

            plotname = paramset_type + suffix
            ax = plots[plotname]
            legend = plotname == 'fixed-mut-inf'
            legend_bbox_to_anchor = (1, 1)
            results_plot(indiv,
                         outcome=outcome,
                         axis=ax,
                         palette=palette,
                         legend=legend,
                         legend_bbox_to_anchor=legend_bbox_to_anchor)

            ax.set_ylabel('')
            ax.set_xlabel('')
            #ax.yaxis.set_major_formatter(ScalarFormatter())
            if 'mut' in plotname:
                ax.set_yscale('log')
                ax.set_yticks([1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
                ax.set_ylim([0.5e-6, 1.5])
            else:
                ax.set_ylim([0, 1.05])


    ## label plot
    plots['fixed-any-inf'].set_title('detectable infections\n'
                                     'per inoculation')
    plots['fixed-mut-inf'].set_title('new variant infections\n'
                                     'per inoculation')
    plots['fixed-mut-inf-per'].set_title('new variant infections\n'
                                     'per detectable infection')
    plots['fixed'].set_ylabel('immediate recall response\n')
    plots['visvar'].set_ylabel('realistic recall response\n')

    plots['all'].set_ylabel('frequency\n\n'
                            ''.format(denom),
                            fontsize='xx-large')
    plots['all'].set_xlabel('bottleneck')

    fig.tight_layout()
    fig.savefig(output_path,
                bbox_inches='tight')

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("\nUSAGE: {} <sensitivity results> "
              "<output path>"
              "\n"
              "".format(sys.argv[0]))
    else:
        main(sys.argv[1],
             sys.argv[2])
