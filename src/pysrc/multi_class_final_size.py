#!/usr/bin/env python3
import numpy as np

def multi_class_final_size(susceptibilities,
                           class_distribution,
                           R0,
                           init_inf=None):
    susceptibilities = np.array(susceptibilities)
    class_distribution = np.array(class_distribution)
    if init_inf is None:
        init_inf = np.zeros_like(class_distribution)
        init_inf[0] = 0.00001
    n_classes = susceptibilities.size
    imm_mat = np.diag(susceptibilities)
    R0_mat = R0 * np.ones_like(imm_mat)
    beta_mat = imm_mat @ R0_mat
    class_distribution = class_distribution / np.sum(class_distribution)
    s_finals = multi_class_S_final(beta_mat,
                                   class_distribution,
                                   init_inf)
    r_finals = class_distribution - s_finals
    return (s_finals, r_finals)

def multi_class_inoculations(susceptibilities,
                             class_distribution,
                             R0,
                             init_inf=None):
    
    s_finals, r_finals = multi_class_final_size(susceptibilities,
                                                class_distribution,
                                                R0,
                                                init_inf)
    Reff = calc_Reff_simple(susceptibilities,
                            class_distribution,
                            R0)
    
    fadeout = min(1 / Reff, 1)
    total_inoculations = (1 - fadeout) * np.sum(R0 * r_finals)
    return (total_inoculations, total_inoculations * np.array(class_distribution))
    
def multi_class_S_final(B_mat, S0s, I0s,
                        recoveries=None,
                        maxiter = 10000):
    if recoveries is None:
        inv_recov = np.identity(S0s.size)
    else:
        inv_recov = np.diag(1 / recoveries)

    def S_iter(X):
        BE = B_mat @ inv_recov
        internal = X - S0s - I0s
        term = BE @ internal
        return S0s * np.exp(term)
    
    X = S0s
    k = 0
    notstop = True
    while(notstop and k < maxiter):
        Xprime = S_iter(X)
        if np.all(np.abs(Xprime - X) < 1e-10):
            notstop = False
        X = Xprime
        k += 1
        
    if notstop:
        raise RuntimeError("Failed to converge in {}"
                           "iterations".format(maxiter))
    return X

def calc_Reff_simple(susceptibilities,
                     class_distribution,
                     R0):
    susceptibilities = np.array(susceptibilities)
    class_distribution = np.array(class_distribution) / np.sum(class_distribution)
    return np.sum(R0 * susceptibilities * class_distribution)

def calc_Reff_full(susceptibilities,
                   class_distribution,
                   R0,
                   recoveries=None):
    class_distribution = np.array(class_distribution)
    class_distribution = class_distribution / np.sum(class_distribution)
    susceptibilities = np.array(susceptibilities)
    imm_mat = np.diag(susceptibilities * class_distribution)
    R0_mat = R0 * np.ones_like(imm_mat)
    beta_mat = imm_mat @ R0_mat
    
    if recoveries is None:
        inv_recov = np.identity(susceptibilities.size)
    else:
        inv_recov = np.diag(1 / recoveries)

    next_gen = beta_mat @ inv_recov
    return np.max(np.linalg.eigvals(next_gen))
