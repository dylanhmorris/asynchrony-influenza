#!/usr/bin/env python3

######################################
# filename: wh_popgen.py
# author: Dylan H. Morris <dhmorris@princeton.edu>
#
# description: functions for analytical
# results and approximations regarding
# the within host competition between
# mutant and wild-type (incl. replication
# selection).
######################################


import parse_parameters as pp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq
from scipy.integrate import quad, dblquad

import point_of_transmission as pt


def freq_mut_select(t, f_0, delta):
    """
    Frequency of a mutant under 
    positive or negative selection
    in the absence of ongoing mutation
    """
    numer = np.exp(delta * t)
    denom = numer + (1 / f_0) - 1
    return numer / denom

def freq_mut_ongoing(t, f_0, delta, mu, growth):
    if not f_0 > 0:
        raise ValueError("initial frequency "
                         "must be greater than "
                         "zero")
    elif not f_0 <= 1:
        raise ValueError("initial frequency "
                         "must be less than "
                         "or equal to one")

    if f_0 == 1:
        return 1

    mg = mu * growth

    if not np.abs(delta) > 0:
        c = 1 / (1 - f_0)
        denom = c + mg * t
        numer = denom - 1
    else:
        c = (delta * f_0 - f_0 * mg + mg) / (f_0 - 1)
        numer = c * np.exp(delta * t) + mg
        denom = numer - delta
        
    return numer / denom


def freq_mut_ongoing_piecewise(
        t_final,
        t_emerge,
        t_M,
        t_peak,
        f_0,
        delta,
        mu,
        growth):
    """
    Handle casework on the timings
    of immunity, ending of mutation
    at t_peak, etc.
    """

    ## determine ordering of t_M, t_peak, t_final, t_emerge

    selection = np.abs(delta) > 0
    
    mutation_ever = t_emerge < t_peak
    mutation_ends = t_peak < t_final

    selection_ever = t_M < t_final and selection
    selection_from_start = t_M <= t_emerge and selection
    
    mutate_and_select_separately = t_M > t_peak and mutation_ever
    mutate_then_mutate_and_select = t_M > t_emerge and mutation_ever

    just_mutate = mutation_ever and not selection_ever

    just_select = selection_ever and not mutation_ever

    do_nothing = not mutation_ever and not selection_ever

    
    if do_nothing:
        result = f_0
        
    elif just_mutate:
        dt_mutate = min(t_peak, t_final) - t_emerge
        
        result = freq_mut_ongoing(
            dt_mutate,
            f_0,
            0,
            mu,
            growth)
        
    elif just_select:
        dt_select = t_final - t_M
        result = freq_mut_select(
            dt_select,
            f_0,
            delta)
        
    elif selection_from_start:
        t_end_mutate = min(t_peak, t_final)
        dt_mutate = t_end_mutate - t_emerge
        dt_select = t_final - t_end_mutate
        
        f_end_mutate = freq_mut_ongoing(
            dt_mutate,
            f_0,
            delta,
            mu,
            growth)
        
        result = freq_mut_select(
            dt_select,
            f_end_mutate,
            delta)
        
    elif mutate_and_select_separately:
        dt_mutate = t_peak - t_emerge
        dt_select = t_final - t_M
        
        f_end_mutate = freq_mut_ongoing(
            dt_mutate,
            f_0,
            0,
            mu,
            growth)
        
        result = freq_mut_select(
            dt_select,
            f_end_mutate,
            delta)
        
    elif mutate_then_mutate_and_select:
        
        t_end_mutate = min(t_peak, t_final)
        dt_just_mutate = t_M - t_emerge
        dt_mutate_and_select = t_end_mutate - t_M
        dt_just_select = t_final - t_end_mutate
        
        f_start_select = freq_mut_ongoing(
            dt_just_mutate,
            f_0,
            0,
            mu,
            growth)
        
        f_end_mutate = freq_mut_ongoing(
            dt_mutate_and_select,
            f_start_select,
            delta,
            mu,
            growth)
        
        result = freq_mut_select(
            dt_just_select,
            f_end_mutate,
            delta)
    else:
        raise ValueError("ordering of events failed!"
                         "Check function and parameters")
        
    return result



def f_m(t_final = None,
        t_M = None,
        delta = None,
        t_emerge = None,
        R0 = None,
        d_v = None,
        f_0 = None,
        mu = None,
        t_peak = None,
        ongoing_mutate = False,
        error_tol = 1e-3):
    
    if t_final is None:
        raise ValueError("Must provide t")
    if t_M is None:
        raise ValueError("Must provide t_M")
    if delta is None:
        raise ValueError("Must provide delta")
    if t_emerge is None:
        raise ValueError("Must provide t_emerge")
    if R0 is None:
        raise ValueError("Must provide R0")
    if d_v is None:
        raise ValueError("Must provide d_v")
    if f_0 is None:
        raise ValueError("Must provide f_0")
    if mu is None:
        raise ValueError("Must provide mu")

    if f_0 == 1:
        ## if at fixation, return 1 
        result = 1
    
    elif t_final < t_emerge:
        ## if not yet emerged, return 0
        result = 0

    elif not ongoing_mutate:
        ## if not accounting for ongoing
        ## mutation or it is 
        t_select = max(0, t_final - max(t_M, t_emerge))
        result = freq_mut_select(
            t_select,
            f_0,
            delta)

    else:
        if t_peak is None:
            raise ValueError("Must provide t_peak "
                             "to handle ongoing mutation")
        ## else need to account
        ## for ongoing mutation
        growth = R0 * d_v

        result = freq_mut_ongoing_piecewise(
            t_final,
            t_emerge,
            t_M,
            t_peak,
            f_0,
            delta,
            mu,
            growth)

    return result


def f_emerge_approx(
        t_emerge,
        t_M,
        k,
        c_w,
        bottleneck,
        R0,
        d_v):

    g0 = R0 * d_v - d_v
    
    if t_emerge < t_M:
        f_0 = (1 / bottleneck) * np.exp(
            -t_emerge * g0)
        
    else:
        f_0 = ((1 / bottleneck) *
               np.exp(-t_M * g0) *
               np.exp(-(t_emerge - t_M) * (g0 - k * c_w)))

    return f_0

def invert_f_m(
        f_target,
        t,
        t_M,
        R0,
        d_v,
        bottleneck,
        k,
        c_w,
        c_m,
        ongoing_mutate = False,
        t_peak = None,
        mu = None):

    if ongoing_mutate:
        if mu is None:
            raise ValueError("Must supply mutation rate "
                             "to invert f_m with ongoing "
                             "mutation")
        if t_peak is None:
            raise ValueError("Must supply time of peak "
                             "to invert f_m with ongoing "
                             "mutation")
    else:
        mu = 0    
    
    def func(t_emerge):
        delta = k * (c_w - c_m)
        freq = f_m(
            t_final = t,
            t_M = t_M,
            delta = delta,
            t_emerge = t_emerge,
            R0 = R0,
            d_v = d_v,
            f_0 = f_emerge_approx(
                t_emerge,
                t_M,
                k,
                c_w,
                bottleneck,
                R0,
                d_v),
            mu = mu,
            t_peak = t_peak,
            ongoing_mutate = ongoing_mutate)

        return freq - f_target
    
    return brentq(func, 0, t)
    

def p_repl_declining(founding_pop,
                     mu,
                     R0,
                     d_v,
                     k,
                     c_w,
                     c_m):
    """
    Calculates an analytical approximate 
    probability of replication selection, given
    no inoculation selection
    """
    Reff_wt = (R0 * d_v) / (d_v + k * c_w)

    ## each founding virion replicates on average
    ## a number of times given by a geometric series
    q = (Reff_wt / (1 - Reff_wt)) 
    Reff_mut = (R0 * d_v) / (d_v + k * c_m)
    if Reff_mut < 1:
        result = 0 
    else:
        p_sse = 1 - (1 / Reff_mut)
        result = q * mu * p_sse * founding_pop

    return result


def t_star_minus(f_target,
                 t,
                 t_M,
                 delta,
                 R0,
                 d_v,
                 bottleneck):


    tau = t - t_M
    g0 = R0 * d_v
    f_ratio = (1 - f_target) / (f_target)
    lnb = np.log(bottleneck)
    logterm = np.log(f_ratio *
                     np.exp(delta * tau) + 1)
    return (logterm - lnb) / (g0 - d_v)

def t_star_plus(f_target,
                t,
                t_M,
                R0,
                d_v,
                bottleneck,
                k,
                c_w,
                c_m):      

    delta = k * (c_w - c_m)
    G0 = R0 * d_v - d_v
    G1 = G0 - k * c_w
    logterm = np.log((1 - f_target) / f_target)
    lnb = np.log(bottleneck)
    numer = logterm - lnb + (delta * t) - (k * c_w * t_M)
    denom = G1 + delta

    return numer / denom

def t_star(f_target,
           t,
           t_M,
           R0,
           d_v,
           bottleneck,
           k,
           c_w,
           c_m):

    delta = k * (c_w - c_m)
    
    t_minus = t_star_minus(
        f_target,
        t,
        t_M,
        delta,
        R0,
        d_v,
        bottleneck)
    
    t_plus = t_star_plus(
        f_target,
        t,
        t_M,
        R0,
        d_v,
        bottleneck,
        k,
        c_w,
        c_m)

    if t_minus < t_M:
        cand_t = max(t_minus, 0)
    else:
        cand_t = max(t_plus, 0)

    return min(t, cand_t)


def p_repl(f_target,
           t,
           t_M,
           R0,
           d_v,
           bottleneck,
           k,
           c_w,
           c_m,
           mu):
    
    if t_M < 0.1 and k * c_w > (R0 - 1) * d_v:
        return p_repl_declining(
            bottleneck,
            mu,
            R0,
            d_v,
            k,
            c_w,
            c_m)
    

    ## get needed emergence
    ## time
    t_star_repl = t_star(
        f_target,
        t,
        t_M,
        R0,
        d_v,
        bottleneck,
        k,
        c_w,
        c_m)
    
    ## get prob of emerging by
    ## that time
            
    result = emergence_time_cdf(t_star_repl,
                                mu,
                                t_M,
                                R0,
                                d_v,
                                k,
                                bottleneck,
                                c_w,
                                c_m)
    return result


    

def emergence_time_cdf(time,
                       mu,
                       t_M,
                       R0,
                       d_v,
                       k,
                       bottleneck,
                       c_w,
                       c_m,
                       verbose = False,
                       return_pdf = False):
    g = R0 * d_v
    b = bottleneck

    def l(t):
        """
        Generation rate per individual as a function
        of time since inoculation
        """
        generation_rate = mu * g
        survival_prob = (g / (g + d_v + k * c_m * (t >= t_M)))
        return generation_rate * survival_prob

    def lamb(t, alpha, pop_0):
        """
        Rate of generation as a function 
        of time and net growth rate (g - d) = alpha
        """
        return pop_0 * np.exp(alpha * t) * l(t)

    def Lamb(t, alpha, pop_0):
        """
        indefinite integral of lamb(t) 
        i.e. the net average events total
        of our variable-rate poisson process
        """
        return (pop_0 * l(t) / alpha) * np.exp(alpha * t)

    if verbose:
        print(l(0))
        print(l(t_M + 1))
        print(Lamb(0, g - d_v, b))
    
    def net_lamb(t):
        """
        definite integral of lamb(x) from 0 to t
        """
        
        if t < t_M:
            return Lamb(t, g - d_v, b) - Lamb(0, g - d_v, b)
        if t >= t_M:
            pop_t_M = b * np.exp((g - d_v) * t_M)
            if verbose:
                print(pop_t_M)
            return (
                Lamb(t_M, g - d_v, b) - Lamb(0, g - d_v, b) +
                Lamb(t - t_M, g - d_v - k * c_w, pop_t_M) -
                Lamb(0, g - d_v - k * c_w, pop_t_M))

    if return_pdf:
        if time < t_M:
            alpha = g - d_v
            gen_rate = b * np.exp(alpha * time) * l(time)
        else:
            alpha_before = g - d_v
            alpha_after = g - d_v - k * c_w
            gen_rate = (b * np.exp(alpha_before * t_M) *
                        np.exp(alpha_after * (time - t_M)) *
                        l(time))
        result = -np.exp(-net_lamb(time)) * (-gen_rate)
        
    else:
        result = 1 - np.exp(-net_lamb(time))
        
    return result


def p_transmit(
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
        t_peak = None,
        p_loss = 0,
        exact = False):
        
    def integrand_func(t_emerge):

        delta = k * (c_w - c_m)

        if t_emerge < t_M:
            f_0 = 1 / (bottleneck * np.exp((R0 - 1) * d_v * t_emerge))
        else:
            V_M = bottleneck * np.exp((R0 - 1) * d_v * t_M)
            f_0 = 1 / (V_M * np.exp((R0 * d_v - d_v - k * c_w) *
                                    (t_emerge - t_M)))
        
        f_transmit = f_m(t_transmit,
                         t_M = t_M,
                         delta = delta,
                         t_emerge = t_emerge,
                         R0 = R0,
                         d_v = d_v,
                         f_0 = f_0,
                         mu = mu,
                         t_peak = t_peak,
                         ongoing_mutate = True)
        return (
            emergence_time_cdf(t_emerge,
                               R0 = R0,
                               d_v = d_v,
                               k = k,
                               t_M = t_M,
                               mu = mu,
                               bottleneck = bottleneck,
                               c_m = c_m,
                               c_w = c_w,
                               return_pdf = True,
                               verbose = False) *
            pt.analytic_p_inocs(bottleneck,
                                f_transmit,
                                p_loss = p_loss,
                                wt_neut_prob = wt_neut_prob,
                                mut_neut_prob = mut_neut_prob,
                                inoculum_size = vb_ratio * bottleneck,
                                exact = exact))

    return quad(integrand_func,
                0,
                t_transmit)[0]


def C_peak(R0, d, C_max):
    g = R0 * d / C_max
    return d / g

def V_repl(r, R0, d, C_max):
    g = R0 * d
    surv = (1 - d/g)
    return r * C_max * (1 - 1/R0)

def mutant_growth_from_mu(Vw, mu, R0, d, C_max=4e8, C=4e8):
    return R0 * d * Vw * mu * C / C_max

def mutant_growth_from_repl(Vm, mu, R0, d, C_max=4e8, C=4e8):
    return R0 * d * Vm * (1 - mu) * C / C_max
