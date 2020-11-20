#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import plotting_style as ps

from point_of_transmission import analytic_p_inocs, poisson_max_sus, neutralize_prob_from_z
from models.susceptibility_models import pick_sus_func

fermi_style = "-."
multiplicative_style = "dotted"
linear_style = "dashed"
titer_style = "solid"
drift_style = "solid"
line_alpha=0.8

def plot_optimal_selector_cluster(
        f_mut = None,
        bottleneck = None,
        R0_wh = None,
        escape = None,
        axis = None,
        legend=True,
        fermi_steepness = None,
        homotypic_titer = None,
        inoculum_size = 25,
        exact = False,
        plot_drift = True,
        inoculum_model = None,
        effective_num_independent = None,
        cross_imm_sigma = None,
        cmaps=[None] * 2,
        darkness = 0.99,
        line_alpha = line_alpha,
        drift_style = drift_style,
        z_homotypic = None,
        verbose = False,
        **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    dists = np.linspace(0, 10, 1000)

    if plot_drift:
        drift_prob = analytic_p_inocs(bottleneck,
                                      f_mut,
                                      p_loss = 0,
                                      mut_neut_prob = 0,
                                      wt_neut_prob = 0,
                                      inoculum_model = inoculum_model,
                                      bottleneck_model = "binomial",
                                      inoculum_size = inoculum_size)
        print("drift emergence prob: {}".format(drift_prob))
        axis.axhline(drift_prob, linestyle=drift_style,
                     color="k", alpha=1)

    models = [
        #"linear",
        "multiplicative",
        #"fermi",
        "titer"
    ]

    styles = [
        # linear_style,
        multiplicative_style,
        # fermi_style,
        titer_style
    ]
    
    for model_name, linestyle, cmap, in zip(
            models,
            styles,
            cmaps):
        print("calculating wild-type and mutant protection "
              "for {} with fermi_steepness {} and "
              "homotypic_titer {}".format(model_name,
                                          fermi_steepness,
                                          homotypic_titer))
        sus_func = pick_sus_func(
            model_name,
            escape,
            fermi_steepness = fermi_steepness,
            homotypic_sus = 1 - z_homotypic)
        
        z_list = [(1 - sus_func(dist),
                   1 - sus_func(dist + 1))
                  for dist in dists]
        
        if effective_num_independent is None:
            effective_num_independent = inoculum_size
            
        wt_neut_probs = [neutralize_prob_from_z(
            z_wt,
            inoculum_size,
            inoculum_model)
                         for z_wt, _ in z_list]
        if cross_imm_sigma is not None:
            mut_neut_probs = [wt_n * cross_imm_sigma
                              for wt_n in wt_neut_probs]
        else:
            mut_neut_probs = [
                neutralize_prob_from_z(
                    z_mut,
                    effective_num_independent,
                    inoculum_model)
                for _, z_mut in z_list]

            
        print("calculating escape probs for {}...".format(model_name))
        probs = [analytic_p_inocs(bottleneck,
                                 f_mut,
                                 p_loss = 0,
                                 mut_neut_prob = mnp,
                                 wt_neut_prob = wnp,
                                 inoculum_size = inoculum_size,
                                 exact = exact,
                                 verbose = False)
                 for mnp, wnp in zip(mut_neut_probs,
                                     wt_neut_probs)]
        if cmap is not None:
            color = cmap(darkness)
        else:
            color = None
        axis.plot(dists, probs,
                  label = model_name,
                  color = color,
                  linestyle = linestyle,
                  alpha = line_alpha,
                  **kwargs)

    if legend:
        axis.legend()
    axis.set_xlabel("cluster distance")
    axis.set_ylabel("emergence probability")


def plot_sus_funcs(axis = None,
                   escape = 0.25,
                   fermi_steepness = 3,
                   fermi_style = fermi_style,
                   multiplicative_style = multiplicative_style,
                   linear_style = linear_style,
                   titer_style = titer_style,
                   cmaps = [None] * 2,
                   max_dark = 0.7,
                   line_alpha = line_alpha,
                   z_homotypic = None,
                   **kwargs):
    if axis is None:
        fig, axis = plt.subplots()

    xs = np.linspace(0, 8, 100)

    for display_name, model_name, model_style, cmap in zip(
            ['multiplicative', 'sigmoid'],
            ['multiplicative', 'titer'],
            [multiplicative_style, titer_style],
            cmaps):
        sus_func = pick_sus_func(model_name,
                                 escape,
                                 fermi_steepness=fermi_steepness,
                                 homotypic_sus = 1 - z_homotypic)
        ys = [sus_func(x) for x in xs]

        if cmap is not None:
            color = cmap(max_dark)
        else:
            color = None
            
        axis.plot(xs, ys, linestyle=model_style,
                  label=display_name,
                  alpha=line_alpha,
                  color=color,
                  **kwargs)

    axis.set_xlabel("cluster distance")
    axis.set_ylabel("susceptibility")
    axis.legend()

def main(fermi_steepness, escape, z_homotypic):
    fig, ax = plt.subplots(1, 2, sharex=True)
    lw=5
    plot_sus_funcs(axis=ax[0],
                   fermi_steepness=fermi_steepness,
                   escape=escape,
                   lw=lw)
    plot_optimal_selector_cluster(axis = ax[1],
                                  fermi_steepness = fermi_steepness,
                                  escape = escape,
                                  z_homotypic = z_homotypic,
                                  lw = lw,
                                  legend = False)
