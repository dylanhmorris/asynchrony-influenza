######################################################
# name: figure_kernel_shift_repl.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: plot how replication selection
# shifts the distribution of fixed mutants
######################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import norm

from wh_popgen import p_transmit
from point_of_transmission import analytic_p_inocs, neutralize_prob_from_z_poisson
from models.susceptibility_models import pick_sus_func



def plot_kernel_shift_repl(
        k = None,
        mu = None,
        d_v = None,
        R0 = None,
        bottleneck = None,
        t_M = None,
        t_transmit = None,
        vb_ratio = None,
        p_loss = 0,
        phenotypes = np.linspace(-1, 1, 250),
        sus_func = None,
        susceptibility_model = 'linear',
        escape = 1,
        wt_phenotype = 0,
        generator_phenotype = None,
        sd_mut = None,
        axis = None,
        z_homotypic = None,
        recipient_phenotype = -99,
        kernel_color = None,
        kernel_linestyle = None,
        kernel_alpha = None,
        dist_color = None,
        dist_alpha = None,
        kernel_lw = None,
        mark_original_phenotype = False,
        wt_phen_lw = 4,
        wt_phen_color = 'black',
        wt_phen_linestyle = 'dashed',
        **kwargs):

    if axis is None:
        fig, axis = plt.subplots()
        
    if generator_phenotype is None:
        generator_phenotype = wt_phenotype

    if sus_func is None:
        sus_func = pick_sus_func(susceptibility_model,
                                 escape = escape,
                                 homotypic_sus = 1 - z_homotypic)
        
    gen_wt_sus = sus_func(abs(wt_phenotype - generator_phenotype))
    recip_wt_sus = sus_func(abs(wt_phenotype - recipient_phenotype))
    c_w = (1 - gen_wt_sus)
    print("generator susceptibility to wild-type: ", gen_wt_sus)
    print("generator wild-type neut rate within host: ", k * c_w)
    recip_wt_neut_prob = neutralize_prob_from_z_poisson(
        (1 - recip_wt_sus) * z_homotypic,
        vb_ratio)

    prob_generate = norm.pdf(phenotypes, wt_phenotype, sd_mut)

    results = []
    for phenotype in phenotypes:
        gen_mut_sus = sus_func(abs(phenotype - generator_phenotype))
        recip_mut_sus = sus_func(abs(phenotype - recipient_phenotype))
        c_m = (1 - gen_mut_sus)
        recip_mut_neut_prob = neutralize_prob_from_z_poisson(
            (1 - recip_mut_sus) * z_homotypic,
            vb_ratio)

        results.append(
            p_transmit(c_w = c_w,
                       c_m = c_m,
                       k = k,
                       mu = mu,
                       d_v = d_v,
                       R0 = R0,
                       bottleneck = bottleneck,
                       t_M = t_M,
                       t_transmit = t_transmit,
                       t_peak = t_transmit,
                       wt_neut_prob = recip_wt_neut_prob,
                       mut_neut_prob = recip_mut_neut_prob,
                       vb_ratio = vb_ratio,
                       p_loss = p_loss))
    prob_survive = np.array(results)
    unnormed_probs = prob_generate * prob_survive
    plot_probs = unnormed_probs / np.sum(unnormed_probs)

    kernel_normed = prob_generate / np.sum(prob_generate)

    ## first plot shifted result, then overlay kernel
    axis.plot(phenotypes,
              plot_probs,
              color = dist_color,
              alpha = dist_alpha,
              **kwargs)

    axis.plot(phenotypes,
              kernel_normed,
              linestyle = kernel_linestyle,
              color = kernel_color,
              alpha = kernel_alpha,
              lw = kernel_lw,
              **kwargs)

    if mark_original_phenotype:
        axis.axvline(wt_phenotype,
                     lw = wt_phen_lw,
                     color = wt_phen_color,
                     linestyle = wt_phen_linestyle)


