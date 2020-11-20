from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.special import lambertw
import numpy as np
import matplotlib.pyplot as plt

def get_q(s0, R0):
    omega = -R0*s0*np.exp(-R0*s0)
    putative_q = float(np.real(lambertw(omega, 0)/(R0*s0) + 1))
    return float(np.max([0, putative_q]))

def get_S_final(s0, R0):
    q = get_q(s0, R0)
    return s0 * (1 - q)

def get_q_bar(s0, c0, R0, c):
    m0 = 1 - s0 - c0
    def temp_func(x):
        val = -R0 * x
        return 1 - x - (
            (s0 * np.exp(val) +
             c0 * np.exp(val * c) + m0))
    return brentq(temp_func, 1e-5, 1-1e-5) 

def get_qs_qc(s0, c0, R0, c):
    q_bar = get_q_bar(s0, c0, R0, c)
    qs = 1 - np.exp(-R0 * q_bar)
    qc = 1 - np.exp(-R0 * q_bar * c)
    return (qs, qc)


def get_mu_bar(s0, c0, R0, c, mu_s, m):
    qs, qc = get_qs_qc(s0, c0, R0, c)
    mu_bar = mu_s * s0 * qs + mu_s * c0 * qc
    return mu_bar


def get_delta_q(s0, R0, N):
    Reff = s0 * R0
    Delta = Reff - 1
    if Delta > 10 * N**(-1/3):
        return get_q(s0, R0)
    elif -Delta > 10 * N**(-1/3):
        qN = 1/(1 - Reff)
        return qN/N
    else:
        return get_mean_final_size(N, s0 * R0)/N


def get_binom_prob(s0, R0, mu, N, alpha = 1):
    Reff = s0 * R0
    q = get_q(s0, R0)
    return 1 - (mu * alpha * (1-s0)) ** (R0 * s0 * q * (1 - (1/Reff)))

    
def get_prob(s0, R0, mu, N, alpha = 1):
    Reff = s0 * R0
    Delta = Reff - 1
    if Delta > 1 * N**(-1/3):
        q = get_q(s0, R0)
        return (1 - np.exp(-(Reff - 1) * q * N * mu * alpha * (1-s0)))
    elif -Delta > 1 * N**(-1/3):
        qN = 1/(1 - Reff)
        return (1 - np.exp(-R0 * qN * mu * alpha * (1-s0)))
    else:
        qN = get_mean_final_size(N, Reff)
        return (1 - np.exp(-R0 * qN * mu * alpha * (1-s0)))


def get_basic_prob(s0, R0, mu, N, alpha = 1):
    Reff = s0 * R0
    q = get_q(s0, R0)
    if q > 0:
        return (1.0 - np.exp(-(Reff - 1) * q * N * mu * alpha * (1-s0)))
    else:
        return 0

def get_small_prob_approx(s0, R0, mu, N, alpha = 1):
    Reff = s0 * R0
    q = get_q(s0, R0)
    if q > 0:
        candidate = ((Reff - 1) * q * N * mu * alpha * (1-s0))
        if candidate < 0.0001:
            return candidate
        else:
            raise ValueError("Probability too large to approximate")
    else:
        return 0

def get_prob_s(S, s0, R0, mu, N, alpha = 1):
    Reff = s0 * R0
    return (1 - (1/Reff)) * (1 - np.exp(- R0 * (s0 - S)
                                        * N * mu * alpha * (1-s0)))


def get_prob_s(S, s0, R0, mu, N, alpha = 1):
    Reff = s0 * R0
    return (1 - np.exp(-((Reff - 1)/s0) * (s0 - S)
                       * N * mu * alpha * (1-s0)))

    
def get_prob_q(s0, R0, mu, N, q):
    return (1 - np.exp(-R0 * s0 * q * N * mu * (1-s0)))


def get_alpha_num(s0, R0, N, m, alpha):
    q = get_q(s0, R0)
    return 1 - np.exp(-R0 *s0 * q *
                      N * m * (1 - alpha) *
                      alpha * (1 - s0))


def get_alpha_num(s0, R0, mu, N, m, alpha):
    q = get_q(s0, R0)
    return 1 - np.exp(-R0 * s0 * q * mu *
                      N * np.exp(-m * alpha) *
                      alpha * (1 - s0))

def womega(R0, S0):
    return lambertw(-R0 * S0 * np.exp(-R0 * S0), 0)

def implicit_s_star(R0, S0):
    w_val = float(np.real(womega(R0, S0)))
    numer = R0 + 1 + np.sqrt(R0**2 + 4 * R0 * w_val**2 +
                             8 * R0 * w_val + 2 * R0 + 1)
    denom = 2 * R0 * (w_val + 2)
    return numer / denom

def get_s_star(R0):
    
    def temp_func(s_val):
        return s_val - implicit_s_star(R0, s_val)
    
    return brentq(temp_func, 1e-5, 1-1e-5)


def get_conditional_s_star(R0):
    
    def temp_func(s_val):
        w_val = float(np.real(womega(R0, s_val)))
        return s_val - (1 / (2 + w_val))
    
    return brentq(temp_func, 1e-5, 1-1e-5)


def s_max(R0):
    def temp_func(s0):
        return -get_prob(s0, R0, 0.000001, 100)

    return brentq(temp_func, 1e-5, 1-1e-5)


def linear_bdk_process(R0, mu, t):
    """
    If there were an infinite and equal supply
    of generators and selectors, what would 
    we get (ignores susceptible depletion)
    """
    b = -(R0 + 1 + mu * R0)
    a = R0
    c = 1
    sqrt_term = np.sqrt(b * b - 4 * a * c)
    v0 = (-b - sqrt_term) / (2 * a)
    v1 = (-b + sqrt_term) / (2 * a)
    v0_expr = np.exp(-R0 * (sqrt_term / a) * t) * (1 - v0)
    numer = (v0 * (v1 - 1)) + (v1 * v0_expr)
    denom = v1 - 1 + v0_expr
    return 1 - (numer / denom)




def get_q_vaccine(R0, theta, I0):
    
    def temp_func(T):
        exp_part = np.exp(-(1 + theta) * R0 * T)
        numer = (1 - theta * T) * (1 + theta) * R0 + theta + (I0 - T) * ((1 + theta)**2) * R0
        denom = (1 + theta) * R0 + theta
        return exp_part - (numer/denom)
    
    return brentq(temp_func, 1e-5, 1-1e-5)

def make_p_s_plot(s0, R0, mu, N, stepsize = 0.001):
    """
    Plot the CDF of susceptible count
    at first mutant (S_mut) from S = s0
    to S_final
    """
    sfinal = s0 - (s0 * get_q(s0, R0))
    s_vals = np.arange(1 - s0, 1 - sfinal, stepsize)
    print(s_vals)
    probs =  [get_prob_s(1 - s_val, s0, R0, mu, N)
              for s_val in s_vals]
    return plt.plot(s_vals, probs)


def get_prob_with_m_depletion(s0, R0, mu, N, alpha = 1):
    Reff = s0 * R0
    Delta = Reff - 1
    if Delta > 1 * N**(-1/3):
        q = get_q(s0, R0)
        return 1 - ((1 - (mu * (1 - np.exp(-(Reff - 1) * q * mu * alpha * (1-s0)))))**N)
    elif -Delta > 1 * N**(-1/3):
        qN = 1/(1 - Reff)
        return mu * (1 - np.exp(-R0 * qN * mu * alpha * (1-s0)))
    else:
        qN = get_mean_final_size(N, Reff)
        return mu * (1 - np.exp(-R0 * qN * mu * alpha * (1-s0)))


def main():
    R0s = np.arange(1.05, 10, 0.05)
    s_maxes = [get_s_star(R0) for R0 in R0s]
    s_conditionals = [get_conditional_s_star(R0) for R0 in R0s]

    plt.plot(R0s, s_maxes)
    plt.plot(R0s, s_conditionals)
    plt.xlabel("R0")
    plt.ylabel("S_opt")
    plt.show()



    N = 10000
    mu = 1e-7
    s0s = np.arange(0, 1.01, 0.01)[::-1]
    R0s = [1.1, 1.3, 1.5, 1.8]
    mylines = []
    fig, ax = plt.subplots(figsize=(8, 8))

    for R0 in R0s:
        xs = 1 - s0s
        approx_probs = [get_basic_prob(s0, R0, mu, N) for s0 in s0s]
        print("finding exact")
        #exact_prob= [np.sum(
        #    get_total_escape_prob(N, int(s0s[5] * N) - 1, R0, mu)[1])]
        #exact_probs = exact_prob * len(s0s)
        mylines += plt.plot(xs, approx_probs, label=str(R0))
        #mylines += plt.plot(s0s, exact_probs, label=str("exact"))
    plt.xlabel("Level of immunity", fontsize=14)
    plt.ylabel("Probability", fontsize=14)
    plt.legend(title="R0 values", handles=mylines, fontsize=14)
    plt.title("Epidemic-level selection (analytical)", fontsize=20)
    fig.savefig("../out/analytic_selection_one.pdf")

    fig, ax = plt.subplots()

    s0s = [0.6, 0.7, 0.8, 0.9]
    mylines = []
    m = 10
    for s0 in s0s:
        alphas = np.arange(0, 1, 0.005)
        alpha_nums = np.array([np.exp(-m*alpha) * get_basic_prob(s0, 5, 1e-10, 500, alpha = alpha) for alpha in alphas])
        alpha_dist = alpha_nums/np.sum(alpha_nums)
        mylines += plt.plot(alphas, alpha_dist, label=str(s0))

    plt.xlabel(r'Jump size $\alpha$')
    plt.ylabel("Frequency")
    plt.title(r'Jump size distribution, $\mu(\alpha) = e^{-10 \alpha}$')
    plt.savefig("../out/intermediate_jump_size_exponential.pdf")
    plt.show()


    s0s = [0.1, 0.5, 0.8]
    mylines = []
    m = 5
    for s0 in s0s:
        alphas = np.arange(0, 1, 0.005)
        alpha_nums = np.array([get_alpha_num(s0, 5, 1e-10, 500, 5, alpha)
                               for alpha in alphas])
        alpha_dist = alpha_nums/np.sum(alpha_nums)
        mylines += plt.plot(alphas, alpha_dist, label=str(s0))

    plt.xlabel("alpha")
    plt.ylabel("p_alpha")
    plt.legend(title="s0 values", handles=mylines)
    plt.show()



    qs = np.arange(0, 1, 0.005)
    probs = np.array([get_prob_q(0.8, 1.05, 1e-5, 100, q) for q in qs])
    plt.plot(qs, probs)
    plt.show()

    myiter = 10000
    periods = 1
    R0s = np.arange(1, 1.5, 0.05)
    results = []
    for R0 in R0s:
        recovery = np.random.exponential(1/periods, myiter)
        for k in range(periods - 1):
            recovery += np.random.exponential(1/periods, myiter)
            
            infection = np.random.exponential(1/R0, myiter)
            
            results.append(np.sum(recovery < infection)/len(infection))



    R0s = np.arange(0, 5, 0.5)
    qs = [(1-(1/R0)) * get_delta_q(1, R0, 1000) for R0 in R0s]
    true_qs = [get_mean_final_size(1000, R0)/1000 for R0 in R0s]
    plt.plot(R0s, qs)
    plt.plot(R0s, true_qs)
    plt.xlabel("R0")
    plt.ylabel("q")
    plt.show()
