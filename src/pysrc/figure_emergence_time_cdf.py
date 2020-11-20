#!/usr/bin/env python3

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import sys
import os

import plotting_style as ps
import parse_parameters as pp
from plotting_functions import setup_multipanel
import seaborn as sns

from wh_popgen import emergence_time_cdf

def get_data(results_dir, bn):
    files = os.listdir(results_dir)
    pattern = "{}.txt".format(bn)

    filenames = [f for f in files if
                 pattern in f]
    if len(filenames) < 1:
        raise ValueError("No file found for pattern {}".format(pattern))
    elif len(filenames) > 1:
        raise ValueError(
            "More than one "
            "file found matching {}".format(pattern))
    
    import_file_path = os.path.join(results_dir, filenames[0])

    return pd.read_table(import_file_path, header=0, delim_whitespace=True)  

def plot_comparison(param_file,
                    model_name,
                    results_dir,
                    bottleneck=5,
                    axis=None,
                    killing=True,
                    label=False):
    """
    Compare simulated model results and 
    prediction from the analytical emergence_time_cdf
    """
    if axis is None:
        fig, axis = plt.subplots()
        
    params = pp.get_params(param_file)

    C_max = pp.get_param("C_max", model_name, params)
    R0 = pp.get_param("R0_wh", model_name, params) 
    d_v = pp.get_param("d_v", model_name, params) 
    mu = pp.get_param("mu", model_name, params) 
    k = pp.get_param("k", model_name, params) 

    data = get_data(results_dir, bottleneck)

    cmap = plt.cm.Reds
    t_Ms = [] #[0, 0.5, 1]
    n_cols = len(t_Ms)
    colors = np.linspace(0.4, 0.9, n_cols)

    if killing:
        bn = 1
    else:
        bn = bottleneck

    times = np.linspace(0, 2, 1000)
    sim_t_M = pp.get_param("t_M", model_name, params)
    c_w = 1
    c_m = 0
    
    probs = [emergence_time_cdf(time,
                                mu,
                                sim_t_M,
                                R0,
                                d_v,
                                k,
                                bottleneck,
                                c_w,
                                c_m)
             for time in times]

    labs = [None, None]
    if label:
        labs = ['simulated', 'analytical']
    emerged = data[data['emergence_time'] < max(times)]['emergence_time']
    sns.distplot(emerged,
                 hist_kws = dict(cumulative = True),
                 kde=False,
                 hist=True,
                 norm_hist=True,
                 bins=np.arange(0, 2, 0.05),
                 ax=axis,
                 label=labs[0])
    axis.plot(times, probs, color="k",
              lw=ps.standard_lineweight,
              label=labs[1])

def setup_figure():
    width = 18.3
    height = (1/2) * width

    lineweight = width / 3
    mpl.rcParams['lines.linewidth'] = lineweight
    mpl.rcParams['axes.formatter.limits'] = (-3, 6)

    # set up multi-panel figure
    fig = plt.figure(figsize=(width, height))
    gs = gridspec.GridSpec(1, 2)


    
    plot_positions = [
        {"name": "fixed",
         "grid_position": np.s_[0, 0],
         "sharex": None,
         "sharey": None},
        
        {"name": "variable",
         "grid_position": np.s_[0, 1],
         "sharex": "fixed",
         "sharey": "fixed"},
                
    ]
    
    plots = setup_multipanel(fig,
                             plot_positions,
                             gridspec=gs,
                             letter_loc=(0.05, 1.07))

    return (fig, plots)

    
def main(fixed_results_dir,
         var_results_dir,
         params_path,
         outpath):

    fixed_model_name = os.path.basename(fixed_results_dir)
    var_model_name = os.path.basename(var_results_dir)

    bottleneck = 1
    
    fig, plots = setup_figure()

    plot_comparison(
        params_path,
        fixed_model_name,
        fixed_results_dir,
        bottleneck,
        axis=plots["fixed"],
        killing=True,
        label=True)

    plots['fixed'].legend()
    
    plot_comparison(
        params_path,
        var_model_name,
        var_results_dir,
        bottleneck,
        axis=plots["variable"],
        killing=True)

    for plot in plots.values():
        plot.set_xlim([-0.05, 2])
        plot.set_ylim([-0.05, 1.05])
        plot.set_xlabel('emergence time', fontsize='xx-large')
        plot.set_ylabel('cumulative probability', fontsize='xx-large')
        plot.label_outer()

    plots['fixed'].set_title('immediate recall response')
    plots['variable'].set_title('delayed recall response')

    fig.tight_layout()
    fig.savefig(outpath)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("\nUSAGE ./figure_emergence_time_cdf.py <fixed immunity data> <variable immunity data> <parameter file> <output_path> \n\n")
    else:
        main(sys.argv[1],
             sys.argv[2],
             sys.argv[3],
             sys.argv[4])
