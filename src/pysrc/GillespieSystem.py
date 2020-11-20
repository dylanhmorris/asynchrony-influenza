import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import root
from scipy.optimize import newton_krylov
import warnings
from sampling import _cached_weighted_choice

class GillespieSystem():

    def __init__(self,
                 propensity_function,
                 initial_state_array,
                 event_matrix=None,
                 initial_time=0,
                 validate=True,
                 event_function=None,
                 clip_lower=-1e100,
                 clip_upper=1e100,
                 cao_epsilon=0.03,
                 negative_tol = -1e-1):

        self.initial_state_array = np.array(initial_state_array)
        self.initial_time = initial_time
        self.propensity_function = propensity_function
        self.clip_lower = clip_lower
        self.clip_upper = clip_upper
        self.cao_epsilon = cao_epsilon
        self.negative_tol = negative_tol
        
        self.ssa_function_dict = {
            "frm": self.first_reaction,
            "direct": self.direct_method
            }

        self.leap_function_dict = {
            "midpoint": self.midpoint_leap,
            "euler": self.euler_leap,
            "grc": self.grc_leap,
            "prc": self.prc_leap,
            "ssi": self.ssi_leap,
            "implicit": self.implicit_leap}

        self.stepsize_function_dict = {
            "prc": self.prc_stepsize,
            "cao": self.cao_stepsize
            }
        
        if event_function is None:
            if event_matrix is None:
                raise RuntimeError(
                    "Either an event matrix (whose ith row "
                    "is the state change given the ith event) "
                    "or a custom event function must be supplied")
            else:
                self.event_matrix = event_matrix
                self.event_function = self.default_event_f
                self.n_events = self.event_matrix.shape[0]
                self.n_vars = self.event_matrix.shape[1]
        else:
            self.event_function = event_function

        if validate:
            # check that user-supplied functions/matrix work
            self.validate_instance_functions()

        self.reset_state()

    def validate_instance_functions(self):
        """
        Runs a couple simple tests 
        to validate that the user-supplied
        propensity_function and event_function/
        event_matrix should work.
        
        By default, called at instatiation
        """
        try:
            test_prop = self.propensity_function(
                self.initial_state_array, 0)
        except:
            raise ValueError(
                "Propensity function must take two arguments: "
                "a state array equal in size to the initial "
                "state array, a time greater than "
                "or equal to zero")
        try:
            test_event = np.zeros(test_prop.size)
            test_array = np.zeros(self.initial_state_array.size)
            for k in range(test_prop.size):
                test_event[k] = 1
                test_array += self.event_function(
                    test_event,
                    self.initial_time)
                test_event[k] = 0
        except:
            raise BaseException(
                "event_function must return an array that can be "
                "added elementwise to the state_array, and must "
                "do so for every event for which the propensity "
                "function provides a rate")
        
    def reset_state(self):
        self.state_array = np.array(self.initial_state_array)
        self.time = self.initial_time
        self.array_ts = [list(self.state_array)]
        self.time_ts = [self.time]
        self.event_ts = []

    def default_event_f(self, event_occurences, time):
        """
        Default event function. Requires that
        a matrix A of update vectors
        (one for each event) is supplied.
        Update vectors should be row-vectors,
        such that a_ij is the change in
        variable j when event i happens.

        Given an array X of event occurences
        (where x_i is the number of times that
        event i has occurred, it returns simply
        the sum over all events i of
        x_i * update_vector i (i.e. the ith
        row of the matrix). This is obatined
        simply by taking the matrix product
        XA
        """
        if type(event_occurences) is int:
            return self.event_matrix[event_occurences]  
        elif np.any(event_occurences < 0):
            raise ValueError("Negative event counts. "
                             "All events must occur "
                             "0 times or more")
        else:
            return event_occurences @ self.event_matrix

    def verify_rates(self, rates):
        """
        Verify that rates are valid numbers
        and at least one is non-zero.

        Returns True if at least one rate is non-zero
        Returns False is all rates are zero
        Raises an exception if any rate is nan
        Raises an exception if any rate is negative
        """
        if np.all(rates == 0):
            return False
        if np.any(np.isnan(rates)):
            raise ValueError("nan rate")
        if np.any(rates < 0):
            raise ValueError("negative rate.\n "
                             "Rates: {0}".format(rates))
        return True

    def get_rates(self, state_array, time):
        return self.propensity_function(state_array, time)
    
    def first_reaction(self, stopping_time):
        rates = self.get_rates(self.state_array, self.time)

        if not self.verify_rates(rates):
            return False

        waiting_times = np.random.exponential(1/rates)

        which_event = waiting_times == waiting_times.min()

        if np.sum(which_event) > 1:
            raise BaseException("Two events have equal min waiting time")

        if self.time + waiting_times[which_event] > stopping_time:
            self.time = stopping_time
            return False

        self.time += waiting_times[which_event][0]  # 0 to keep time
        # as a float rather than array

        self.event_ts.append(np.argmax(which_event))

        self.state_array += self.event_function(which_event, self.time)
        return True

    def direct_method(self, stopping_time):
        """
        Take a step according to the direct method
        of Gillespie
        """
        rates = self.get_rates(self.state_array, self.time)
        if not self.verify_rates(rates):
            return False

        total_rate = np.sum(rates)
    
        waiting_time = np.random.exponential(1 / total_rate)

        if self.time + waiting_time > stopping_time:
            self.time = stopping_time
            return False

        self.time += waiting_time
        which_event = _cached_weighted_choice(rates, total_rate)
        self.event_ts.append(which_event)
        
        self.state_array += self.event_function(which_event, self.time)
        
        return True

    def calc_approx_midpoint(self, state_array, rates, timestep):
        """
        Calculates an approximate midpoint
        when there is a static, time- and state-array
        invariant event_matrix
        """
        rates_of_change = rates @ self.event_matrix
        midpoints = state_array + 0.5 * timestep * rates_of_change
        return midpoints

    def deterministic_df(self, state_array, time):
        """
        Return the deterministic rates of change
        for the system at a given state array
        """
        return self.get_rates(state_array, time) @ self.event_matrix

    def midpoint_leap(self, timestep):
        """
        Performs a midpoint tau-leap step
        """
        rates = self.get_rates(self.state_array, self.time)
        if not self.verify_rates(rates):
            return False

        approx_midpoint = self.calc_approx_midpoint(
            self.state_array,
            rates,
            timestep)

        zero_vars = approx_midpoint < 0
        if np.any(zero_vars):
            warnings.warn("Negative midpoint...setting variable to zero but "
                          "simulation behavior may be erratic")
            self.state_array[zero_vars] = 0
            approx_midpoint[zero_vars] = 0
        event_rates = timestep * self.get_rates(
            approx_midpoint,
            self.time + 0.5 * timestep)

        event_counts = np.random.poisson(event_rates)

        self.time += timestep
        self.state_array += self.event_function(event_counts, self.time)

        # truncate leaps to below 0, if this is asked for
        if np.any(self.state_array < self.negative_tol):
            print(self.state_array)
            raise RuntimeError("Leap to below negative tolerance!"
                               "Try reducing tau or using an "
                               "adaptive step size")
        self.state_array = np.clip(self.state_array,
                                   self.clip_lower,
                                   self.clip_upper)
        return True

    def eta(self, state_array, t):
        result = np.zeros(self.n_events ** 2).reshape((self.n_events, -1))
        rates = self.get_rates(state_array, t)
        for k in range(self.n_events):
            a_js_k = self.get_rates(state_array + self.event_matrix[k], t)
            result[:,k] = a_js_k - rates
        return result 

    def prc_leap(self, timestep):
        eta = self.eta(self.state_array, self.time)
        rates = self.get_rates(self.state_array, self.time)
        means = rates * timestep + (
            0.5 * timestep * timestep * eta @ rates)
        zero_vars = means < 0
        if np.any(zero_vars):
            warnings.warn("Negative mean occurences...setting to zero but "
                          "simulation behavior may be erratic\n\n"
                          "Try setting tau smaller\n")
            means[zero_vars] = 0
        r_star = np.random.poisson(means)
        self.time += timestep
        self.state_array += r_star @ self.event_matrix
        self.state_array = np.clip(self.state_array,
                                   self.clip_lower,
                                   self.clip_upper)

    def grc_leap(self, timestep):
        warnings.warn("GRC leap not currently implemented properly")
        eta = self.eta(self.state_array, self.time)
        rates = self.get_rates(self.state_array, self.time)
        r = np.random.poisson(rates * timestep)
        E_r = 0.5 * timestep * eta @ r
        Var_r = 0.5 * timestep * timestep * eta @ rates
        # fix negative variance terms
        for j in range(self.n_events):
            k = 0
            while(Var_r[j] < 0):
                if eta[j, k] < 0:
                    E_r[j] += 0.5 * timestep * (
                        (rates[k]/rates[j]) * r[j] - timestep * rates[k])
                    
                    Var_r[j] -= timestep * timestep * rates[k] * eta[j, k]
                k += 1
        r_corr = np.random.normal(E_r, Var_r)

        self.time += timestep
        self.state_array += (r + r_corr) @ self.event_matrix
        self.state_array = np.clip(self.state_array,
                                   self.clip_lower,
                                   self.clip_upper)

    def euler_leap(self, timestep):
        """
        Take a step using Euler tau leaping
        bc why not?
        """
        rates = self.get_rates(self.state_array, self.time)
        r = np.random.poisson(rates * timestep)
        self.time += timestep
        self.state_array += r @ self.event_matrix
        self.state_array = np.clip(self.state_array,
                                   self.clip_lower,
                                   self.clip_upper)

    def ssi_leap(self, timestep):
        z = self.state_array
        def imp_func(y):
            y_rates = self.get_rates(y, self.time) 
            return z + (timestep * y_rates @ self.event_matrix) - y
        y_guess = self.get_rates(z, self.time) * timestep @ self.event_matrix
        y = root(imp_func, y_guess, method="lm")
        rates = self.get_rates(y.x, self.time + timestep)
        self.state_array += np.random.poisson(rates * timestep) @ self.event_matrix
        self.time += timestep

    def implicit_leap(self, timestep):
        x_old = self.state_array
        rates = self.get_rates(x_old, self.time)
        random_term = np.random.poisson(rates * timestep) - rates * timestep
        x_guess = x_old + (rates * timestep) @ self.event_matrix

        def implicit_equation(x_new):
            occurences = random_term + self.get_rates(x_new, self.time + timestep) * timestep
            delta_x = occurences @ self.event_matrix
            return x_old + delta_x - x_new

        soln = newton_krylov(implicit_equation, x_guess)

        if np.any(soln < self.negative_tol):
            print(soln)
            raise RuntimeError("Negative solution "
                               "beyond rounding tolerance! "
                               "Try decreasing tau or using "
                               "an adaptive algorithm")

        self.state_array = np.clip(soln,
                                   self.clip_lower,
                                   self.clip_upper)
        self.time += timestep

    def cao_stepsize(self):
        rates = self.get_rates(self.state_array, self.time)
        mean_i = rates @ self.event_matrix
        var_i = rates @ (np.square(self.event_matrix))
        return 0

    def prc_stepsize(self, max_tau):
        rates = self.get_rates(self.state_array, self.time)
        eta = self.eta(self.state_array, self.time)
        eta_cond = rates / (10 * ((eta @ rates) + 1e-10))
        candidate_taus = np.abs(eta_cond[eta_cond < 0])
        if len(candidate_taus) > 0:
            tau = min(np.min(candidate_taus), max_tau)
        else:
            tau = max_tau
        return tau
    
    def simulate_timeseries(self, max_iter, stopping_time,
                            method = "direct"):
        iter = 0
        not_stop = True
        # flag to check that system not at equilibrium
        # and self.time is less than stopping_time

        step = self.ssa_function_dict[method]
        
        while (iter < max_iter and not_stop):
            not_stop = step(stopping_time)
            # copy array to save snapshot
            self.array_ts.append(np.array(self.state_array))
            self.time_ts.append(self.time)
            iter += 1

            
    def simulate_adaptive_leap_timeseries(
            self, stopping_time, max_tau, max_iter=5000000,
            leap_type="prc", stepsize_algorithm="prc"):
        iter = 0
        not_stop = True
        # flag to check that system not at equilibrium
        # and self.time is less than stopping_time

        leap_function = self.leap_function_dict[leap_type]

        stepsize_function = self.stepsize_function_dict[stepsize_algorithm]
        
        while (iter < max_iter and self.time < stopping_time):
            timestep = stepsize_function(max_tau)
            leap_function(timestep)
            # copy array to save snapshot
            self.array_ts.append(np.array(self.state_array))
            self.time_ts.append(self.time)
            iter += 1


    def simulate_leap_timeseries(self, n_iter, timestep,
                                 leap_type = "midpoint"):

        leap_function = self.leap_function_dict[leap_type]

        for iter in range(int(n_iter)):
            leap_function(timestep)
            self.array_ts.append(list(self.state_array))
            self.time_ts.append(self.time)

    def integrate_timeseries(self, stopping_time, n_steps):
        times = np.linspace(self.time, stopping_time, n_steps)
        
        timeseries = odeint(self.deterministic_df,
                            self.state_array,
                            times)
        self.array_ts += list(timeseries)
        self.time_ts += list(times)
        self.time = stopping_time
        self.state_array = self.array_ts[-1]
        
    def plot_ts(self, variables, logged=False):
        if not logged:
            for col in variables:
                plt.plot(self.time_ts,
                         [item[col] for item in self.array_ts])
        else:
            for col in variables:
                plt.plot(self.time_ts,
                         [np.log10(item[col] + 1) for
                          item in self.array_ts])
                       
    def __repr__(self):
        return "GillespieSystem with state {0}".format(
            self.state_array)


def ones_event_matrix(n_vars):
    """
    Returns a standard event matrix
    for a system in which each variable
    is affected by two events -- an increment
    event and a decrement event -- each of
    which affects no other variables.

    The ordering of the matrix is 0++, 0--, 1++, 
    1--, etc, for the 0th, 1st, 2nd, variable etc
    """
    result = np.zeros(2 * n_vars * n_vars).reshape((2 * n_vars, -1))
    for k in range(n_vars):
        result[2 * k, k] = 1
        result[2 * k + 1, k] = -1
    return result 
