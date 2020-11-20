#!/usr/bin/env python3

# Calculate the final size distribution for a standard stochastic SIR
# epidemic
#
# Adapted from the algorithm and MATLAB code of Black and Ross,
# "Computation of epidemic final size distributions", Journal of
# Theoretical Biology, 2015. https://doi.org/10.1016/j.jtbi.2014.11.029

from scipy.misc import comb 
import numpy as np
from scipy.stats import binom
import matplotlib.pyplot as plt

def psi(N, k):
    return comb(N + k + 1, N, exact=True)

def get_final_size_probs(N, R0, I_initial = 1):
    beta = R0 / (N - 1)
    gamma = 1
    q = np.zeros(N + 1)
    q[I_initial] = 1

    for Z2 in range(0, N + 1):
        for Z1 in range(Z2 + 1, N):
            p1 = 1 / ( 1 + gamma / (beta * ( N - Z1)))
            q[Z1 + 1] = q[Z1 + 1] + q[Z1] * p1
            q[Z1] = q[Z1] * (1 - p1)

    return q


def get_mean_final_size(N, R0, I_initial = 1):
    q = get_final_size_probs(N, R0, I_initial)
    return np.sum(q * np.arange(0, N + 1))


def get_total_escape_prob(N, S0, R0, mu, I0 = 1):
    Qpop = S0 + I0
    M0 = N - Qpop
    invbeta = (N - 1)/R0
    q = np.zeros(N + 1)
    q[I0] = 1
    p_e = np.zeros(N + 1)

    for Z2 in range(0, Qpop + 1):
        for Z1 in range(Z2 + 1, Qpop + 1):
            S = Qpop - Z1
            z1_plus = S
            z2_plus =  invbeta
            escape = mu * M0 
            total_rate = z1_plus + z2_plus + escape
            if(z1_plus > 0):
                q[Z1 + 1] = q[Z1 + 1] + q[Z1] * z1_plus / total_rate
            p_e[Z1] += q[Z1] * escape / total_rate
            q[Z1] = q[Z1] * z2_plus / total_rate

    return (q, p_e)


def p_e_drift(N, S0, R0, mu, bottleneck, I0 = 1):
    mu_s = 1 - binom.pmf(0, bottleneck, mu)
    mu_d =  1 - binom.cdf(bottleneck / 2,
                          bottleneck,
                          mu)
    print(mu_s, mu_d)
    Qpop = S0 + I0
    M0 = N - Qpop
    invbeta = (N - 1)/R0
    q = np.zeros(N + 1)
    q[I0] = 1
    p_e = np.zeros(N + 1)

    for Z2 in range(0, Qpop + 1):
        for Z1 in range(Z2 + 1, Qpop + 1):
            S = Qpop - Z1
            z1_plus = S
            z2_plus =  invbeta
            escape = mu_s * M0 + mu_d * S
            total_rate = z1_plus + z2_plus + escape
            if(z1_plus > 0):
                q[Z1 + 1] = q[Z1 + 1] + q[Z1] * z1_plus / total_rate
            p_e[Z1] += q[Z1] * escape / total_rate
            q[Z1] = q[Z1] * z2_plus / total_rate

    return (q, p_e)



def get_total_escape_prob_loss(N, S0, R0, mu, I0 = 1):
    Qpop = S0 + I0
    M0 = N - Qpop
    beta = R0 / (N - 1)
    gamma = 1
    q = np.zeros(N + 1)
    q[I0] = 1
    p_e = np.zeros(N + 1)

    for Z2 in range(0, Qpop + 1):
        for Z1 in range(Z2 + 1, Qpop + 1):
            S = Qpop - Z1
            z1_plus = (beta * S)
            z2_plus =  (gamma)
            escape = beta * mu * M0 * np.max([0, (1-(N/(R0*(M0 + S))))])
            total_rate = z1_plus + z2_plus + escape
            if(z1_plus > 0):
                q[Z1 + 1] = q[Z1 + 1] + q[Z1] * z1_plus / total_rate
            p_e[Z1] += q[Z1] * escape / total_rate
            q[Z1] = q[Z1] * z2_plus / total_rate

    return (q, p_e)


def escape_prob_m_depletion(N, S0, R0, mu, I0 = 1):
    Qpop = S0 + I0
    M0 = N - Qpop
    beta = R0 / (N - 1)
    gamma = 1
    q = np.zeros(N + 1)
    q[I0] = 1
    p_e = np.zeros(N + 1)

    for Z2 in range(0, Qpop + 1):
        for Z1 in range(Z2 + 1, Qpop + 1):
            S = Qpop - Z1
            z1_plus = (beta * S)
            z2_plus =  (gamma)
            escape = beta * mu * np.max([(M0 - Z1), 0])
            total_rate = z1_plus + z2_plus + escape
            if(z1_plus > 0):
                q[Z1 + 1] = q[Z1 + 1] + q[Z1] * z1_plus / total_rate
            p_e[Z1] += q[Z1] * escape / total_rate
            q[Z1] = q[Z1] * z2_plus / total_rate

    return (q, p_e)


def plot_p_e_bottleneck(N, S0, R0, mus,
                        bottlenecks=np.arange(1, 15)):
    lines = []
    for mu in mus:
        probs = [np.log10(np.sum(p_e_drift(N, S0, R0, mu, bottleneck)[1]))
                 for bottleneck in bottlenecks]
        lines += plt.plot(bottlenecks, probs, label=str(mu))
    plt.xlabel("Bottleneck size")
    plt.ylabel("Log P_mutant")
    plt.legend(title="mu", handles=lines)

def plot_p_e_S0(N, S0s, R0, mu, I0 = 1):
    plt.plot(S0s, [np.sum(get_total_escape_prob(N, S0, R0, mu, I0)[1])
                   for S0 in S0s])
