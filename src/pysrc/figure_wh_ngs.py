#!/usr/bin/env python3

##############################################
# name: figure_wh_ngs.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: summary of within-host next-gen
# sequencing studies
##############################################

import matplotlib.pyplot as plt
import plotting_style as ps
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
from plotting_functions import subdivided_hist

def subst_type(data,
               min_freq = 0):
    freq = data["freq_var"]
    if_ridge = data["is_antigenic_ridge"]
    if_site = data["is_antigenic_site"]
    if_subst = data["is_aa_substitution"]
    
    if if_ridge and freq > min_freq:
        return "ridge"
    elif if_site and freq > min_freq:
        return "site"
    elif if_subst and freq > min_freq:
        return "subst"
    else:
        return "no subst"

def plot_ngs_hist(data,
                  bins = np.arange(0, 1, 0.05),
                  site_color =  ps.antigenic_site_color,
                  ridge_vax_color = ps.ridge_vax_color,
                  ridge_no_vax_color = ps.antigenic_ridge_color,
                  subst_color = ps.generic_substitution_color,
                  vax_column = "vaccination_matched",
                  axis = None,
                  legend = False,
                  min_freq = 0):
    
    if axis is None:
        fig, axis = plt.subplots()

    non_ant_subst = ((data['is_aa_substitution'])
                     & (~data['is_antigenic_site']))
    subst_hist, _= np.histogram(
        data[non_ant_subst]['freq_var'],
        bins = bins)

    site_only = ((data['is_antigenic_site'])
                 & (~data['is_antigenic_ridge']))
    site_hist, _= np.histogram(
        data[site_only]['freq_var'],
        bins = bins)

    ridge_vax = data['is_antigenic_ridge'] & data[vax_column]
    ridge_vax_hist, _= np.histogram(
        data[ridge_vax]['freq_var'],
        bins = bins)

    ridge_no_vax = (data['is_antigenic_ridge'] &
                    ~data[vax_column])
    ridge_no_vax_hist, _= np.histogram(
        data[ridge_no_vax]['freq_var'],
        bins = bins)

    countvecs = np.array([
        ridge_vax_hist,
        ridge_no_vax_hist,
        site_hist,
        subst_hist])

    colorvec = np.array([
        ridge_vax_color,
        ridge_no_vax_color,
        site_color,
        subst_color])
    
    labelvec = np.array([
        'ridge,\nvax.',
        'ridge',
        'site',
        'non-ant.'])
    

    to_use = np.array(
        [np.sum(item) > 0 for item in countvecs])
    
    subdivided_hist(countvecs[to_use],
                    bins = bins,
                    normed = False,
                    colors = colorvec[to_use],
                    labels = labelvec[to_use],
                    axis = axis,
                    legend = legend)
    if legend:
        handles, labels = axis.get_legend_handles_labels()
        axis.legend(handles = handles,
                    labels = labels)

    axis.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
    axis.set_xlim([0, 0.5])
    
def plot_wh_ngs(data,
                axis = None,
                min_freq = 0,
                vax_column = "vaccination_class",
                site_color = ps.antigenic_site_color,
                ridge_color = ps.antigenic_ridge_color,
                subst_color = ps.generic_substitution_color,
                nothing_color = ps.no_substitution_color,
                nothing_alpha = 0.75,
                bar_width = 0.95,
                legend = False,
                **kwargs):
    
    if axis is None:
        fig, axis = plt.subplots()
    
    data['subst_type'] = data.apply(
        subst_type,
        axis = 1)

    cats_in_order = [
        "subst",
        "site",
        "ridge",
        "no subst"]

    ## must be in same order
    ## as categories, to make
    ## pandas barplot of crosstab
    ## work as expected
    barcolors =  [
        subst_color,
        site_color,
        ridge_color,        
        nothing_color]


    data['subst_type_cat'] = pd.Categorical(
        data['subst_type'],
        categories = cats_in_order,
        ordered = True)

    zeros = data[~data["is_aa_substitution"]]
    positives = data[data["is_aa_substitution"]]
    
    zct = pd.crosstab(
        zeros[vax_column],
        zeros["subst_type_cat"],
        dropna = False)

    pct = pd.crosstab(
        positives[vax_column],
        positives["subst_type_cat"],
        dropna = False)

    zeros_colors = [nothing_color]
    
    zct.plot.bar(
        ax = axis,
        stacked = True,
        color = zeros_colors,
        legend = False,
        alpha = nothing_alpha,
        width = bar_width,
        rot = 0,
        **kwargs)

    if not pct.empty:
        print(pct)
        pct.plot.bar(
            ax = axis,
            stacked = True,
            color = barcolors,
            legend = False,
            alpha = 1,
            width = bar_width,
            rot = 0,
            **kwargs)

    if legend:
        handles = [
            mpatches.Patch(
                facecolor = subst_color,
                alpha = 1,
                label = 'non-ant.',
                **kwargs),
            mpatches.Patch(
                facecolor = site_color,
                alpha = 1,
                label = 'site',
                **kwargs),
            mpatches.Patch(
                facecolor = ridge_color,
                alpha = 1,
                label = 'ridge',
                **kwargs),
            mpatches.Patch(
                facecolor = nothing_color,
                alpha = 1,
                label = 'none',
                **kwargs)]

        axis.legend(handles = handles,
                    labels = [h.get_label() for h in handles],
                    loc = "upper left",
                    fancybox = True,
                    frameon = True,
                    handlelength = 0.5,
                    fontsize = 'x-small')

    axis.set_xlabel('vaccination status')
    axis.set_xticklabels(['match',
                          'mismatch',
                          'none'],
                         fontsize = 'small')


def make_presentation_plot():
    """
    makes a version of the figure
    for presentations
    """
    empirical_data_file = "../../dat/cleaned/wh_data.csv"
    empirical_data = pd.read_csv(empirical_data_file)
    width = 15
    height = 7
    fig, axes = plt.subplots(1, 2,
                             figsize=(width, height))
    plot_wh_ngs(empirical_data,
                axis = axes[0],
                legend = True)
    plot_ngs_hist(empirical_data,
                  axis = axes[1],
                  min_freq = 0,
                  legend = True)
    axes[0].set_ylabel("number of infections")
    axes[0].set_xlabel("vaccination status")
    axes[1].set_xlabel("variant frequency")
    axes[0].set_title('observed antigenic substitutions')
    axes[1].set_title('variant within-host frequencies')
    fig.savefig("../../cons/wh_ngs_presentation.pdf")
