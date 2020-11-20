#!/usr/bin/env python3

###################################################
# filename: figure_prepl_analytic_vs_sims.py
# author: Dylan H. Morris <dhmorris@princeton.edu>
# 
# description: generates SI figure
# showing analytial vs simulated probabilities
# of replication selection
###################################################

import sys
import os

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

import plotting_style as ps
from wh_popgen import p_repl

mpl.rcParams['lines.markeredgewidth'] = 1

def calc_quantities(x):
    p_repl_one_percent = np.mean((x.end_mutant_freq >= 0.01))
    p_repl_consensus = np.mean((x.end_mutant_freq >= 0.5))    
    
    return (
        pd.Series(
            [
                p_repl_one_percent,
                p_repl_consensus
            ],
        
            index=['one-percent', 'consensus'])
        )


def get_data(path_to_data,
             pattern=".txt"):
    """
    Import data separated across
    a bunch of flat files
    """
    files = os.listdir(path_to_data)
    df = pd.DataFrame()
    for f in files:
        if re.search(pattern, f):
            filepath = os.path.join(path_to_data, f)
            print("Reading ", filepath)
            temp_dat = pd.read_table(filepath,
                                     delim_whitespace=True)
            df = pd.concat([df, temp_dat], axis=0)
    return df

def main(path_to_data=None,
         output_path=None):

    if path_to_data is None:
        path_to_data = '../../out/replication_selection_results/prepl.txt'
        
    if output_path is None:
        output_path = (
            '../../ms/supp/figures/figure-prepl-analytic-vs-sims.pdf')
    
    cmap = plt.cm.Purples

    print('reading in data...')
    dat = get_data(path_to_data)

    dat['end_mutant_freq'].fillna(value = 0,
                                  inplace=True)

    print('calculating simulated quantities from data...')
    results = dat.groupby(
        ['k',
         't_M']).apply(calc_quantities)

    results = pd.melt(results.reset_index(),
                      id_vars=['k',
                               't_M'],
                      var_name='outcome',
                      value_name='prob')

    print('setting up figure...')

    height = 18.3
    width = (2 / 3) * height

    fig, axes = plt.subplots(ncols=2,
                             nrows=pd.unique(dat.k).size,
                             sharex=True,
                             sharey=True,
                             figsize=(width, height))

    t_Ms = np.linspace(0, 3, 250)
    ks = np.sort(pd.unique(dat.k))
    print(ks)
    print(axes)

    for ax_row, k in enumerate(ks):
        for ax_col, outcome, f_target in zip(
                [0, 1],
                ['one-percent', 'consensus'],
                [0.01, 0.5]):
        
            df = results[
                (results.outcome == outcome) &
                (results.k == k)]
        
            color = ps.repl_color

            t_final = 3
            bottleneck = 1
            R0 = 5
            d_v = 4
            c_w = 1
            c_m = 0
            mu = .33e-5
            
            analytical = [
                p_repl(f_target,
                       t_final,
                       t_M,
                       R0,
                       d_v,
                       bottleneck,
                       k,
                       c_w,
                       c_m,
                       mu)
                for t_M in t_Ms]
        
            axes[ax_row, ax_col].plot(
                t_Ms,
                analytical,
                linestyle = 'dashed',
                lw=width * 0.4,
                color = color,
                alpha = 0.8)

            axes[ax_row, ax_col].plot(
                df.t_M,
                df.prob,
                marker='o',
                color = color,
                markersize=width,
                linestyle = "",
                markeredgecolor='k',
                markeredgewidth=width * 0.1,
                label = outcome)

            axes[ax_row, ax_col].set_yscale('log')

            axes[ax_row, ax_col].set_ylabel('probability for k = {}'.format(k))
            axes[ax_row, ax_col].set_xlabel('$t_M$')
            axes[ax_row, ax_col].label_outer()        

    # style and save plot
    axes[0, 0].set_title('one percent')
    axes[0, 1].set_title('consensus')
    fig.tight_layout()
    fig.savefig(output_path)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("\nUSAGE: {} <path to data> <output path>"
              "\n"
              "".format(sys.argv[0]))

    else:
        main(path_to_data = sys.argv[1],
             output_path = sys.argv[2])
