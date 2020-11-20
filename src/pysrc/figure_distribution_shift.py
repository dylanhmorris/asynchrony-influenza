#!/usr/bin/env python3

# plots how the mucosal barrier
# shifts the distribution of virions

import matplotlib.pyplot as plt
import numpy as np
import plotting_style as ps
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

import analyze_selection_comparison as asc
from point_of_transmission import get_pre_filter, get_post_filter, get_pre_filter_full_joint
from plotting_functions import subdivided_hist                                   
def inoc_class(wt_virions, mut_virions):
    if wt_virions > 0 and mut_virions > 0:
        return 'both'
    elif wt_virions > 0:
        return 'wild-type'
    elif mut_virions > 0:
        return 'mutant'
    else:
        return 'neither'

vec_inoc_class = np.vectorize(inoc_class)


def squarify(M, val, min_dim=0):
    """
    With thanks to user
    yannis on StackOverflow 
    https://stackoverflow.com/questions/10871220/making-a-matrix-square-and-padding-it-with-desired-value-in-numpy
    """
    (a, b)=M.shape
    if a < min_dim and b < min_dim:
        padding = ((0, min_dim - a), (0, min_dim - b))
    elif a > b:
        padding = ((0, 0), (0, a - b))
    else:
        padding=((0, b - a), (0, 0))
    return np.pad(M,
                  padding,
                  mode='constant',
                  constant_values=val)


def plot_sim_inocula(dataframe,
                     inoculum_sizes=None,
                     axis=None,
                     wt_inoc_col="Vw_inoc",
                     mut_inoc_col="Vm_inoc",
                     group_col = "n_encounters_iga",
                     f_mut_threshold=0,
                     min_dim=0,
                     colors={'neither': ps.no_inf_color,
                             'wild-type': ps.wt_color,
                             'mutant': ps.mut_color,
                             'both': ps.mixed_inoc_color}):

    if axis is None:
        fig, axis = plt.subplots()
        
    mut_founding = dataframe.peak_mutant_freq > f_mut_threshold
    subset = dataframe[mut_founding]
    if inoculum_sizes is not None:
        subset = subset[subset.n_encounters_iga.isin(inoculum_sizes)]

    subset['inoc_class'] = vec_inoc_class(subset[wt_inoc_col],
                                          subset[mut_inoc_col])

    print(subset[['inoc_class', 'Vw_iga', 'Vm_iga']])
    
    groups = (
        (subset.groupby(group_col)['inoc_class']\
         .value_counts()) /
        (subset.groupby(group_col)['inoc_class'].size()))
    print(groups)

    label=True
    for v in inoculum_sizes:
        bottom = 0
        for name in ['neither', 'wild-type', 'mutant', 'both']:
            val = groups[v].get(name, 0)
            if label:
                lab = name
            axis.bar(str(v),
                     val,
                     bottom=bottom,
                     edgecolor='k',
                     linewidth=1.25,
                     color=colors.get(name),
                     label=lab)
            bottom += val
        label = False
    

def plot_marginal_dist(full_joint,
                       axis=None,
                       logged=True,
                       which="mutant",
                       **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    sum_axis_dict = {
        "mutant": 1,
        "wt": 0
        }

    sum_axis = sum_axis_dict[which]
    
    naive_bottleneck = full_joint.shape[0] - 1
    n_mutants = np.arange(naive_bottleneck + 1)
    p_mutants = np.sum(full_joint, axis=sum_axis)
    axis.plot(n_mutants,
              p_mutants,
              "o",
              **kwargs)
    if logged:
        axis.set_yscale("log")
        axis.set_ylim(ymax=2)


def plot_percent_mutant(full_joint,
                        axis=None,
                        logged=True,
                       **kwargs):
    if axis is None:
        fig, axis = plt.subplots()
    
    naive_bottleneck = full_joint.shape[0] - 1
    n_mutants = np.arange(naive_bottleneck + 1)

    
    p_mutants = np.sum(full_joint, axis=1)
    axis.plot(n_mutants,
              p_mutants,
              "o",
              **kwargs)
    if logged:
        axis.set_yscale("log")
        axis.set_ylim(ymax=2)

    
def inoculum_heatmap(full_joint,
                     wt_min=0,
                     mut_min=0,
                     wt_max=None,
                     mut_max=None,
                     min_prob=1e-7,
                     logged=True,
                     cmap=plt.cm.Purples,
                     axis=None):
    if axis is None:
        fig, axis = plt.subplots()

    naive_bottleneck = full_joint.shape[0] - 1
    if wt_max is None:
        wt_max = naive_bottleneck
    if mut_max is None:
        mut_max = naive_bottleneck
        
    extent = [wt_min - 0.5, wt_max + 0.5,
              mut_max + 0.5, mut_min - 0.5]
    subset = full_joint[mut_min:mut_max + 1,
                        wt_min:wt_max + 1]

    if logged:
        heatmap_norm = colors.LogNorm(
            vmin=min_prob,
            vmax=1)
    else:
        heatmap_norm=colors.Normalize(
            vmin=0.,vmax=1.)

    heatmap = axis.imshow(subset,
                          cmap=cmap,
                          norm=heatmap_norm,
                          alpha=1,
                          extent=extent,
                          interpolation='nearest',
                          origin='upper',
                          aspect=1)
    axis.invert_yaxis()
    axis.set_xlabel("wild-type virions")
    axis.set_ylabel("mutant virions")
    #axis.set_xticks(np.arange(full_joint.shape[0] + 1))
    #axis.set_yticks(np.arange(full_joint.shape[1] + 1))

    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(heatmap, cax=cax)


def get_percentages(full_joint):
    percentages = np.zeros_like(full_joint)
    for n_mut in range(naive_bottleneck):
        for n_wt in range(naive_bottleneck):
            if n_mut + n_wt > 0:
                percentages[n_mut, n_wt] = n_mut / (n_mut + n_wt)
            else:
                percentages[n_mut, n_wt] = 0
    return percentages

def mutant_hist(full_joint,
                **kwargs):
    percentages = get_percentages(full_joint)
    return np.histogram(percentages,
                        weights=full_joint,
                        **kwargs)


def plot_filtered_inocula(inoculum_sizes,
                          f_mutants,
                          axis=None,
                          kappa_w = 0,
                          kappa_m = 0,
                          colors=[ps.no_inf_color,
                                  ps.wt_color,
                                  ps.mut_color,
                                  ps.mixed_inoc_color],
                            **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    inoculum_sizes = np.array(inoculum_sizes)
    f_mutants = np.array(f_mutants)
        
    both_vals = (
        (1 - np.exp(-inoculum_sizes * (1 - kappa_w) * (1 - f_mutants))) *
        (1 - np.exp(-inoculum_sizes * (1 - kappa_m) * f_mutants))
    )
    
    wt_vals = (1 - np.exp(-inoculum_sizes *
                          (1 - kappa_w) *
                          (1 - f_mutants)) -
               both_vals)
    mut_vals = (1 - np.exp(-inoculum_sizes *
                           (1 - kappa_m) *
                           f_mutants) -
                both_vals)
    
    neither_vals = np.exp(-inoculum_sizes * (
        (1 - kappa_m) * f_mutants + (1 - kappa_w) * (1 - f_mutants)))
    
    all_vals = [neither_vals,
                wt_vals,
                mut_vals,
                both_vals]
    val_names = ['no infection',
                 'wild-type',
                 'mutant',
                 'both']


    bottom = np.zeros_like(all_vals[0])
    
    for vals, color, val_name in zip(all_vals,
                                     colors,
                                     val_names):
        axis.bar([str(item) for item in inoculum_sizes],
                 vals,
                 bottom=bottom,
                 color=color,
                 linewidth=1.25,
                 edgecolor='k',
                 label=val_name)
        bottom += vals
