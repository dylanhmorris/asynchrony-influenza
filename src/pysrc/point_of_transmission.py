#!/usr/bin/env python3

######################################
# filename: point_of_transmission.py
# author: Dylan H. Morris <dhmorris@princeton.edu>
#
# description: functions for analytical
# results and approximations regarding
# the point of transmission / inoculation
# selection
######################################

from models.susceptibility_models import pick_sus_func
from scipy.stats import binom, poisson, hypergeom
import numpy as np
import warnings
from multi_class_final_size import multi_class_inoculations, multi_class_final_size


def p_cib(w,
          b):
    uncorrected = (
        b *
        (1 - np.exp(-w)) /
        w)
            
    correction_term = np.sum([
        (1 - (b / (x + 1))) *
        poisson.pmf(x, w)
        for x in range(b)])
    return uncorrected + correction_term


def analytic_p_inocs(bottleneck,
                     f_mut,
                     p_loss = 0,
                     wt_neut_prob = None,
                     mut_neut_prob = None,
                     exact = True,
                     inoculum_model = "poisson",
                     bottleneck_model = "binomial",
                     inoculum_size = None,
                     get_type_probs = False,
                     verbose = False):
    """
    Calculates an analytical approximate 
    probability of inoculation selection
    with inoculum prior to neutralization either
    equal to naive bottleneck (default)
    or larger, with a bottlenecking afterward
    """
    if inoculum_size is None:
        inoculum_size = bottleneck

    mut_post_filter_mean = inoculum_size * f_mut * (1 - mut_neut_prob)

    p_mutant_survives_mucosa = 1 - np.exp(-mut_post_filter_mean)
    
    wt_post_filter_mean = inoculum_size * (1 - f_mut) * (1 - wt_neut_prob)

    if exact:
        max_mut = poisson.ppf(1-1e-10, mut_post_filter_mean)
        mut_integration_range = np.arange(0, max_mut + 1)
        
        max_wt = poisson.ppf(1-1e-10, wt_post_filter_mean)
        wt_integration_range = np.arange(0, max_wt + 1)

        wt_survival_probs = poisson.pmf(wt_integration_range, wt_post_filter_mean)
        mut_survival_probs = poisson.pmf(mut_integration_range, mut_post_filter_mean)
        post_bn_dist = post_bottleneck_dist(
            wt_survival_probs,
            mut_survival_probs,
            bottleneck)
        p_mutant_survives_bottleneck = np.sum(post_bn_dist[1:,:])
        p_mutant_survives_mucosa = np.sum(mut_survival_probs[1:])

        if get_type_probs:
            return post_bn_dist
        else:
            pass
        
    else:
        if wt_neut_prob < 1:
            uncorrected = (
                bottleneck *
                (1 - np.exp(-wt_post_filter_mean)) /
                wt_post_filter_mean)
            
            correction_term = np.sum([
                (1 - (bottleneck / (x + 1))) *
                poisson.pmf(x, wt_post_filter_mean)
                for x in range(bottleneck)])
            p_mutant_survives_competitive_bottleneck = (
                uncorrected + correction_term)
            # for numerical stability
            if verbose:
                print("uncorrected prob: ", uncorrected)
                print("correction_term: ", correction_term)
                print("p survives competitive bottleneck: ",
                      p_mutant_survives_competitive_bottleneck)

        else:
            p_mutant_survives_competitive_bottleneck = 1
            
        p_mutant_survives_bottleneck = (
            p_mutant_survives_mucosa *
            p_mutant_survives_competitive_bottleneck)
    
    # p(mutant survives) is sum of all entries in the dist
    # with non-zero mutant.
    if verbose:
        print("wt neutralization prob: ", wt_neut_prob)
        print("mut neutralization prob: ", mut_neut_prob)
        print("wt mean post filter: ", wt_post_filter_mean)
        print("prob mutant survives mucosa: ", p_mutant_survives_mucosa)
        print("prob mutant survives bottleneck: ", p_mutant_survives_bottleneck)
    result = (
        p_mutant_survives_bottleneck * (1 - p_loss))

    return result 



def minimal_p_loss(f_mut, R0, bottleneck):
    mean_mut_inoculum = max(f_mut * bottleneck, 1)
    return (1 / R0) ** mean_mut_inoculum

def npzp(z, v):
    return 1 + np.log(z * (1 - np.exp(-v)) + np.exp(-v)) / v

def weighted_minimal_p_inocs(
        bottleneck,
        f_mut,
        host_frequencies=[1],
        wt_neut_probs=[1],
        mut_neut_probs=[0],
        inoculum_size=None,
        exact=False):
    host_frequencies = np.array(host_frequencies) / np.sum(host_frequencies)
    if np.any(host_frequencies < 0):
        raise ValueError("Cannot have negative frequencies")
    weighted_ps = [host_freq * analytic_p_inocs(bottleneck,
                                       f_mut,
                                       p_loss = 0,
                                       wt_neut_prob=wt_neut,
                                       mut_neut_prob=mut_neut,
                                       inoculum_size=inoculum_size,
                                       exact=exact)
                   for host_freq, wt_neut, mut_neut in
                   zip(host_frequencies,
                       wt_neut_probs,
                       mut_neut_probs)]
    return np.sum(weighted_ps)
    
def get_pre_filter(f_mutant, naive_bottleneck):
    n_mutants = np.arange(0, naive_bottleneck + 1)
    p_n_mutants = binom.pmf(n_mutants, naive_bottleneck, f_mutant)
    return (n_mutants, p_n_mutants)

def get_pre_filter_full_joint(f_mutant, naive_bottleneck):
    n_mutants, p_n_mutants = get_pre_filter(f_mutant, naive_bottleneck)
    return np.diagflat(p_n_mutants[::-1])[::-1]

            
def get_post_filter(f_mutant,
                    inoculum_size,
                    z_wt,
                    z_mut,
                    inoculum_model="poisson"):
    mut_pre, p_pre = get_pre_filter(f_mutant, inoculum_size)

    wt_neut_prob = neutralize_prob_from_z(
        z_wt, inoculum_size, inoculum_model)
    mut_neut_prob = neutralize_prob_from_z(
        z_mut, inoculum_size, inoculum_model)

    joint = np.zeros((inoculum_size + 1) ** 2).reshape(
        (inoculum_size + 1, -1))

    for n_m_pre, p_m_pre in zip(mut_pre, p_pre):
        n_w_pre = inoculum_size - n_m_pre
        for n_m_post in range(n_m_pre + 1):
            for n_w_post in range(n_w_pre + 1):
                cond_p_mw_post = (
                    binom.pmf(n_m_post,
                              n_m_pre,
                              1 - mut_neut_prob) *
                    binom.pmf(n_w_post,
                              n_w_pre,
                              1 - wt_neut_prob))
                joint[n_m_post, n_w_post] += (
                    p_m_pre * cond_p_mw_post)
    return joint


def bottleneck_joint_dist(joint_dist, bottleneck):
    if type(joint_dist) is not np.ndarray:
        raise ValueError("joint_dist must by numpy array")
    
    outcome = np.zeros((bottleneck + 1) ** 2).reshape(
        (bottleneck + 1, -1))
    joint_size = joint_dist.shape[0]
        
    for n_m_pre in range(joint_size):
        for n_w_pre in range(joint_size):
            if n_m_pre + n_w_pre < bottleneck:
                outcome[n_m_pre, n_w_pre] += joint_dist[n_m_pre, n_w_pre]
            else:
                p_pre = joint_dist[n_m_pre, n_w_pre]
                m_freq = n_m_pre / (n_w_pre + n_m_pre)
                for n_m_post in range(bottleneck + 1):
                    cond_p_n_m_post = binom.pmf(n_m_post,
                                                bottleneck,
                                                m_freq)
                    n_w_post = bottleneck - n_m_post
                    outcome[n_m_post, n_w_post] += (
                        p_pre * cond_p_n_m_post)
    return outcome


def z_poisson(mean_inoculum, neutralization_prob):
    """
    Returns the probability that an infection is fully
    neutralized if inocula are Poisson with
    mean mean_inoculum and virions are 
    independently neutralized with probability
    neutralization_prob.
    """
    return np.exp(-mean_inoculum * (1 - neutralization_prob))

def neutralize_prob_from_z_poisson(
        desired_z,
        mean_inoculum,
        conditional_on_inoculated = True):
    """
    Returns an implied mucosal
    neutralization probability for a given
    probability of full neutralization desired_z
    and a given mean poisson distributed
    inoculum size mean_inoculum
    """
    
    if conditional_on_inoculated:
        prob_nonzero_naive = 1 - np.exp(-mean_inoculum)
        desired_sus = prob_nonzero_naive * (1 - desired_z)
        desired_z = 1 - desired_sus

    result = np.zeros_like(desired_z)
    result[desired_z != 0] = 1 + np.log(desired_z) / mean_inoculum
    if np.any(result < 0):
        warnings.warn("needed neutralization prob was"
                      "negative, returning 0")
        print("desired z: ", desired_z)
        print("mean_inoculum v:", mean_inoculum)
        print(result)
    result = np.maximum(0, result)
    return result 

def neutralize_prob_from_z_fixed(
        desired_z,
        inoculum,
        conditional_on_inoculated = True):
    """
    Returns an implied mucosal
    neutralization probability for a given
    probability of full neutralization desired_z
    and a given fixed inoculum size
    inoculum
    """

    if inoculum < 1:
        raise ValueError("inoculum size must be at least 1")
    return desired_z ** (1 / inoculum)

def neutralize_prob_from_z(
        desired_z,
        inoculum_size,
        inoculum_model,
        conditional_on_inoculated = True):

    models = {
        "fixed": neutralize_prob_from_z_fixed,
        "poisson": neutralize_prob_from_z_poisson
    }

    model = models.get(inoculum_model)
    if model is None:
        raise KeyError("No inoculum model {}. Available models:"
                       "{}".format(inoculum_model,
                                   [key for key in models.keys()]))
    return model(
        desired_z,
        inoculum_size,
        conditional_on_inoculated = conditional_on_inoculated)

def poisson_max_sus(inoculum_size):
    """
    In the poisson model, 
    there's a maximum susceptibility
    given inocula because some inocula
    would result in 0 virions P(0) 
    
    this mathematical convenience 
    saves us from having to 
    zero-truncate the poisson
    for other analyses
    """
    return 1 - np.exp(-inoculum_size)


def post_bottleneck_dist(
        wt_count_at_bn_probs,
        mut_count_at_bn_probs,
        bottleneck):
    results = np.zeros((bottleneck + 1)**2).reshape(
        (bottleneck + 1, bottleneck + 1))
    for n_muts, p_n_muts in enumerate(mut_count_at_bn_probs):
        for n_wts, p_n_wts in enumerate(wt_count_at_bn_probs):
            total = n_muts + n_wts
            if total > bottleneck:
                for n_muts_final in range(bottleneck+ 1):
                    results[n_muts_final, bottleneck - n_muts_final] += (
                        p_n_muts * p_n_wts * 
                        hypergeom.pmf(n_muts_final, total, n_muts, bottleneck))
            else:
                results[n_muts, n_wts] += p_n_muts * p_n_wts
    return results


def inocs_by_pct_immune(
        f_mut = None,
        final_bottleneck = None,
        mucus_bottleneck = None,
        escape = None,
        fermi_steepness = None,
        mut_neut_prob = None,
        wt_neut_prob = None,
        x_min = 0,
        x_max = 1,
        fineness = 1000,
        inoculum_model = None,
        bottleneck_model = None):

    wt_pcts = np.linspace(x_min, x_max, fineness)
    if x_max > 1:
        raise ValueError("Cannot have more than 100% immune")

    
    per_inoc_inocs_prob = [
        wt_pct * analytic_p_inocs(
            final_bottleneck,
            f_mut,
            p_loss = 0,
            wt_neut_prob = wt_neut_prob,
            mut_neut_prob = mut_neut_prob,
            inoculum_size = mucus_bottleneck,
            inoculum_model = inoculum_model,
            bottleneck_model = bottleneck_model,
            exact = False) +
        (1 - wt_pct) * analytic_p_inocs(
            final_bottleneck,
            f_mut,
            p_loss = 0,
            wt_neut_prob = 0,
            mut_neut_prob = 0,
            inoculum_size = mucus_bottleneck,
            inoculum_model = inoculum_model,
            bottleneck_model = bottleneck_model,
            exact = False)
        for wt_pct in wt_pcts]
    
    return (wt_pcts, per_inoc_inocs_prob)


def get_inocs_from_immunity(
        f_mut = None,
        final_bottleneck = None,
        mucus_bottleneck = None,
        escape = None,
        fermi_steepness = None,
        distribution = [],
        x_min = 0,
        x_max = 1,
        susceptibility_model = None,
        cross_imm_sigma = None,
        homotypic_sus = None,
        fineness = 1000,
        inoculum_model = "poisson"):

    sus_func = pick_sus_func(susceptibility_model,
                             escape,
                             fermi_steepness = fermi_steepness,
                             homotypic_sus = homotypic_sus)
    
    z_wts = get_z_wts(distribution, sus_func)

    wt_neut_probs = np.array([
        neutralize_prob_from_z(
            z_wt,
            mucus_bottleneck,
            inoculum_model)
        for z_wt in z_wts])
    
    z_muts = get_z_muts(distribution, sus_func)

    if cross_imm_sigma is not None:
        mut_neut_probs = cross_imm_sigma * wt_neut_probs
    else:
        mut_neut_probs = np.array([
            neutralize_prob_from_z(
                z_mut,
                mucus_bottleneck,
                inoculum_model)
            for z_mut in z_muts])

    print(wt_neut_probs)
    print(weighted_minimal_p_inocs(
        final_bottleneck,
        f_mut,
        wt_neut_probs=[wt_neut_probs[0]],
        mut_neut_probs=[mut_neut_probs[0]],
        inoculum_size=mucus_bottleneck))
    wt_pcts = np.linspace(x_min, x_max, fineness)
    if x_max > 1:
        raise ValueError("Cannot have more than 100% immune")
    per_inoc_inocs_prob = [
        weighted_minimal_p_inocs(final_bottleneck,
                                 f_mut,
                                 host_frequencies=get_freqs(wt_pct,
                                                            distribution),
                                 wt_neut_probs=wt_neut_probs,
                                 mut_neut_probs=mut_neut_probs,
                                 inoculum_size=mucus_bottleneck)
        for wt_pct in wt_pcts]

    return (wt_pcts, per_inoc_inocs_prob)

def get_z_wts(distribution, sus_func):
    return np.array([1 - sus_func(0)] +
                    [1 - sus_func(x) for x, _ in distribution] +
                    [1 - sus_func(99)])

def get_z_muts(distribution, sus_func):
    return np.array([1 - sus_func(1)] +
                    [1 - sus_func(x + 1) for x, _ in distribution] +
                    [1 - sus_func(99)])



def get_freqs(wt_pct, distribution):
    freqs = np.array([0] + [y for _, y in distribution] + [0])
    remaining = 1 - np.sum(freqs)
    freqs[-1] = remaining
    if remaining < 0:
        raise ValueError("distribution of known phenotypes "
                         "cannot sum to more than 1")
    elif any([freq < 0 for freq in freqs]):
        raise ValueError("distribution of known phenotypes "
                         "must only include non-negative "
                         "frequencies 0 <= f <= 1")

    n_freqs = len(freqs)
    if n_freqs > 2:
        new_freqs = np.array([x for x in freqs])
        new_freqs[0] = wt_pct
        to_go = np.sum(new_freqs) - 1
        k = 1
        while to_go > 0 and k < n_freqs + 1:
            change = min(to_go, freqs[-k])
            new_freqs[-k] -= change
            to_go = np.sum(new_freqs) - 1
            k += 1
    else:
        new_freqs = np.array([wt_pct, 1 - wt_pct])
    return new_freqs


def epidemic_p_inocs(
        wt_pct,
        f_mut = None,
        bottleneck = None,
        v = None,
        R0_population = None,
        escape = None,
        fermi_steepness = None,
        distribution = [],
        susceptibility_model = "linear",
        inoculum_model = "poisson",
        cross_imm_sigma = None,
        z_wt = None,
        exact = False,
        N = 1):
    
    sus_func = pick_sus_func(susceptibility_model,
                             escape,
                             fermi_steepness,
                             homotypic_sus = 1 - z_wt)

    z_wts = get_z_wts(distribution, sus_func)
    wt_susses = 1.0 - z_wts # tell numpy this is a float!
    freqs = get_freqs(wt_pct, distribution)

    tot_inoc, class_inoc = multi_class_inoculations(
        wt_susses,
        freqs,
        R0_population)

    s_finals, r_finals = multi_class_final_size(wt_susses,
                                                freqs,
                                                R0_population)
    imm_final = np.sum(r_finals)
    s_finals = list(s_finals) + [imm_final]
        
    log_per_cap_none = [
        (class_inoc *
         np.log(1 - analytic_p_inocs(bottleneck,
                                     f_mut,
                                     p_loss = 0,
                                     wt_neut_prob = neutralize_prob_from_z(
                                         z_wt,
                                         v,
                                         inoculum_model,
                                         conditional_on_inoculated = True),
                                     mut_neut_prob = neutralize_prob_from_z(
                                         z_wt,
                                         v,
                                         inoculum_model) * cross_imm_sigma,
                                     inoculum_size = v,
                                     exact=exact)))
        for z_wt, class_inoc in zip(z_wts, class_inoc)]

    return 1 - np.exp(np.sum(log_per_cap_none * N))
