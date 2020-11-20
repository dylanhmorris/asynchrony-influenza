#!/usr/bin/env python3

import numpy as np
from scipy.optimize import fsolve

default_titer_alpha = 2.844
default_titer_beta = 1.299

def fermi_susceptibility(distance, steepness=1, half=2, max_sus=1):
    const = 1 + np.exp(-steepness * half)
    return 1 - (const / (1 + np.exp(steepness * (distance - half))))


def get_implied_ln_titer(
        alpha,
        beta,
        sus):
    ln_titer = alpha + np.log(1/sus - 1) / beta
    return ln_titer

def get_fold_drop_per_cluster(
        escape,
        ln_homotypic_titer,
        alpha=default_titer_alpha,
        beta=default_titer_beta):
    sus_of_1 = escape
    ln_needed_titer = get_implied_ln_titer(alpha,
                                           beta,
                                           sus_of_1)
    needed_titer = np.exp(ln_needed_titer)
    print("Implied one cluster titer: ", needed_titer)
    ln_fold_drop = ln_homotypic_titer - ln_needed_titer
    result = np.exp(ln_fold_drop)
    print("Implied fold drop per cluster: ", result)
    return result 
    
def dist_to_sus_via_titer(distance,
                          fold_drop_per_cluster=2**5,
                          mean_homotypic_titer=300,
                          **kwargs):
    predicted_titer = mean_homotypic_titer/fold_drop_per_cluster**distance
    return susceptibility_curve(predicted_titer)


def susceptibility_curve(titer,
                         alpha=default_titer_alpha,
                         beta=default_titer_beta):
    """
    Implements the estimated protection
    curve by HI titer from 
    Coudeville et al 2010
    https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-10-18
    but as a susceptibility curve
    """
    return (1 / (1 + np.exp(beta * (np.log(titer) - alpha))))


def inv_fermi_sus(sus, steepness=1, half=2):
    const = 1 + np.exp(-steepness * half)
    exponent = np.log((1 - const - sus) / (sus - 1))
    return (exponent / steepness) + half

def get_fermi_half(sus_of_1, steepness):
    def func_solve(x):
        return fermi_susceptibility(1, steepness, x) - sus_of_1
    return fsolve(func_solve, 1)[0]
    
def linear_susceptibility(distance, escape, homotypic_sus):
    return min(1, homotypic_sus + distance * (escape - homotypic_sus))

def multiplicative_susceptibility(distance,
                                  escape,
                                  homotypic_sus):
    return 1 - (1 - homotypic_sus) * ((1 - escape)/(1 - homotypic_sus))**distance

def pick_sus_func(susceptibility_model,
                  escape,
                  fermi_steepness = None,
                  ln_homotypic_titer = None,
                  homotypic_sus = 0,
                  alpha = default_titer_alpha,
                  beta = default_titer_beta):
    
    if susceptibility_model == "fermi":
        fermi_half = get_fermi_half(escape, fermi_steepness)
        def sus_func(dist):
            return fermi_susceptibility(
                dist, fermi_steepness, fermi_half)
        
    elif susceptibility_model == "titer":
        if ln_homotypic_titer is None:
            ln_homotypic_titer = get_implied_ln_titer(alpha,
                                                      beta,
                                                      homotypic_sus)
        fdpc = get_fold_drop_per_cluster(
            escape,
            ln_homotypic_titer = ln_homotypic_titer)
        mean_ht = np.exp(ln_homotypic_titer)
        
        def sus_func(dist):
            return dist_to_sus_via_titer(
                dist,
                fold_drop_per_cluster = fdpc,
                mean_homotypic_titer = mean_ht)
    else:
        sus_funcs = {
            "multiplicative": multiplicative_susceptibility,
            "linear": linear_susceptibility,
        }
        def sus_func(dist):
            return sus_funcs[susceptibility_model](dist,
                                                   escape,
                                                   homotypic_sus)
        
    return sus_func
