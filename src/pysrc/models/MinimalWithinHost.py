# MinimalWithinHost.py defines a stochastic model
# of within-host flu dynamics that adapts
# a simple target cell limited model to exert
# a non-constant degree of adaptive
# immune pressure

from GillespieSystem import GillespieSystem
import numpy as np

class MinimalWithinHost(GillespieSystem):
    """
    Implements a simple target-cell-virion  
    kinetics model with a potentially 
    non-constant adaptive response and 
    with antigenic mutation
    """

    def __init__(self,
                 beta,
                 r_w,
                 r_m,
                 mu,
                 d_v,
                 k,
                 t_M,
                 Cmax,
                 initial_state_array,
                 adaptive_response_type="binary",
                 c=0,
                 clip_lower=0,
                 p_M = None,
                 f = 1,
                 p_C = 0,
                 E_w = 0,
                 E_m = 0,
                 t_N = 9999
                 ):
        self.beta = beta
        self.r_w = r_w
        self.r_m = r_m
        self.f = f
        self.mu = mu
        self.d_v = d_v
        self.k = k
        self.t_M = t_M
        self.Cmax = Cmax
        self.p_M = p_M
        self.p_C = p_C
        self.c = c
        self.adaptive_response_type = adaptive_response_type
        self.E_w = E_w
        self.E_m = E_m
        self.t_N = t_N
        
        events = np.array(
        [[ 1,  0,  0 ],
         [-1,  0,  0 ],
         [ 0,  1,  0 ],
         [ 0, -1,  0 ],
         [ 0,  0,  1 ],
         [ 0,  0, -1 ]])

        self.m_funcs = {"binary": self.binary_adaptive_response,
                        "logistic": self.logistic_adaptive_response}
        
        super(MinimalWithinHost, self).__init__(
            None,
            initial_state_array,
            events,
            validate=False,
            clip_lower=clip_lower)

    def get_rates(self, state_vector, t):

        C = state_vector[0]
        Vw = state_vector[1]
        Vm = state_vector[2]
        Mw_eff = max(self.E_w, self.c * self.E_m)
        Mm_eff = max(self.E_m, self.c * self.E_w)

        C_plus = self.p_C * (C + 1) * (1 - (C / self.Cmax)) * (C <= self.Cmax)
        C_minus =  self.f * self.beta * C * (Vw + Vm)

        Vw_plus = self.r_w * self.beta * C * (Vw + self.mu * Vm)
        Vw_minus = (self.beta * C * Vw) + (self.d_v + self.k * max(Mw_eff * self.M(t), self.N(t))) * Vw
        
        Vm_plus =  self.r_m * self.beta * C * (Vm + self.mu * Vw) 
        Vm_minus = (self.beta * C * Vm) + (self.d_v + self.k * max(Mm_eff * self.M(t), self.N(t))) * Vm

        rate_array = np.clip(np.array([C_plus, C_minus,
                                       Vw_plus, Vw_minus,
                                       Vm_plus, Vm_minus]),
                             0, 1e100)
        
        return rate_array

    def M(self, t):
        return self.m_funcs[self.adaptive_response_type](t)
    
    def N(self, t):
        return t > self.t_N

    def binary_adaptive_response(self, t):
        return t > self.t_M

    def logistic_adaptive_response(self, t):
        x = self.p_M * (t - self.t_M)
        return np.exp(x) / (1 + np.exp(x))


def calc_beta(R0, Cmax, d_v, r):
    beta = R0 * d_v / (r * Cmax)
    return beta
