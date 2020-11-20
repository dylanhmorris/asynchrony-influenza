#!/usr/bin/env python3

######################################################
# name: figure_relative_contribution.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: plot relative role of replication
# and inoculuation selection in increasing
# mutant frequency
######################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import norm
import sys

import wh_popgen as whp
import plotting_style as ps
import plotting_functions as pf
import point_of_transmission as pt


def get_probs(
        k = None,
        c_w = None,
        c_m = None,
        mu = None,
        d_v = None,
        R0 = None,
        bottleneck = None,
        t_M = None,
        t_transmit = None,
        wt_neut_prob = None,
        mut_neut_prob = None,
        vb_ratio = None,
        t_peak = None):

    p_drift = whp.p_transmit(
        k = 0,
        c_w = 0,
        c_m = 0,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = bottleneck,
        t_M = t_M,
        t_transmit = t_transmit,
        wt_neut_prob = 0,
        mut_neut_prob = 0,
        vb_ratio = vb_ratio,
        t_peak = t_peak,
        p_loss = 0)

    p_repl_only = whp.p_transmit(
        k = k,
        c_w = c_w,
        c_m = c_m,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = bottleneck,
        t_M = t_M,
        t_transmit = t_transmit,
        wt_neut_prob = 0,
        mut_neut_prob = 0,
        vb_ratio = vb_ratio,
        t_peak = t_peak,
        p_loss = 0)

    p_inoc_only = whp.p_transmit(
        k = 0,
        c_w = 0,
        c_m = 0,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = bottleneck,
        t_M = t_M,
        t_transmit = t_transmit,
        wt_neut_prob = wt_neut_prob,
        mut_neut_prob = mut_neut_prob,
        vb_ratio = vb_ratio,
        t_peak = t_peak,
        p_loss = 0)

    p_all = whp.p_transmit(
        k = k,
        c_w = c_w,
        c_m = c_m,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = bottleneck,
        t_M = t_M,
        t_transmit = t_transmit,
        wt_neut_prob = wt_neut_prob,
        mut_neut_prob = mut_neut_prob,
        vb_ratio = vb_ratio,
        t_peak = t_peak,
        p_loss = 0)

    return [p_drift,
            p_inoc_only,
            p_repl_only,
            p_all]

    
def main(outpath):

    delta_taus = np.linspace(0, 3, 10)
    
    k = 20
    mu = 0.33e-5
    d_v = 4
    R0 = 5
    bottleneck = 1
    t_M = 2
    t_transmit = 3
    wt_neut_prob = 1
    z_m = 0.75
    vb_ratios = [10, 100, 1000]
    mut_neut_probs = [0, 0.5, 0.99]
    
    t_peak = 2
    tau = 0.5

    ## figure setup
    n_rows = len(vb_ratios)
    n_cols = len(mut_neut_probs)
    width = 15
    height = 15
    figsize = (width, height)
    
    fig, axes = plt.subplots(n_rows,
                             n_cols,
                             figsize = figsize,
                             sharex = True,
                             sharey = True)
    colors = [
        "k",
        ps.inocs_color,
        ps.repl_color,
        ps.both_color]

    linestyles = [
        "dashed",
        "solid",
        "dotted",
        "-."]
        
    labels = [
        "drift",
        "inoc",
        "repl",
        "both"]

    for ind, vb_ratio in enumerate(vb_ratios):
        for mut_neut_prob, ax in zip(mut_neut_probs, axes[ind,:]):
            results = np.array([
                get_probs(
                    k = k,
                    c_w = 1,
                    c_m = 1 - ((delta_tau / tau) / k),
                    mu = mu,
                    d_v = d_v,
                    R0 = R0,
                    bottleneck = bottleneck,
                    t_M = t_M,
                    t_transmit = t_M + tau,
                    wt_neut_prob = wt_neut_prob,
                    mut_neut_prob = mut_neut_prob,
                    vb_ratio = vb_ratio,
                    t_peak = t_peak)
                for delta_tau in delta_taus])

        
            for x, color, ls, label in zip(range(4),
                                           colors,
                                           linestyles,
                                           labels):
                ax.plot(delta_taus,
                        results[:, x],
                        color = color,
                        linestyle = ls,
                        label = label)
                if ind < 1:
                    ax.set_title("$\kappa_m = {}$"
                                 "".format(mut_neut_prob))
        
                ax.set_yscale('log')
                ax.set_xlim(left = 0, right = 3)

    axes[0, 0].legend()

    for ax in axes[:, 2]:
        ax.yaxis.set_label_position("right")
    axes[0, 2].set_ylabel("$v = 10$",
                          labelpad = 25,
                          rotation = 270,
                          fontsize = "x-large")
    axes[1, 2].set_ylabel("$v = 100$",
                          labelpad = 25,
                          rotation = 270,
                          fontsize = "x-large")
    axes[2, 2].set_ylabel("$v = 1000$",
                          labelpad = 25,
                          rotation = 270,
                          fontsize = "x-large")

    bound = pf.add_bounding_subplot(fig)

    bound.set_ylabel("probability of new variant "
                     "transmission $p_\mathrm{nv}$",
                     labelpad = 25,
                     fontsize = "xx-large")
    bound.set_xlabel("transmitting host selection $\delta \\tau$",
                     labelpad = 25,
                     fontsize = "xx-large")

    fig.tight_layout()

    fig.savefig(outpath)



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("USAGE: ./{} <outpath>".format(sys.argv[0]))
    else:
        main(sys.argv[1])
