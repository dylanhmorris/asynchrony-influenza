#!/usr/bin/env python3

########################################################
# filename: plot_wh_timecourse.py
# author: Dylan Morris <dhmorris@princeton.edu
# description: provides helper function for plotting
# how simulated timecourses of infection according
# to the within-host model (both absolute and virion
# counts and variant frequencies). Can also super-
# impose analytical predictions
########################################################


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import analyze_selection_comparison as asc
import parse_parameters as pp
from plotting_functions import plot_df_timecourse, plot_df_freqs, scale_color_brightness
from wh_popgen import f_m
import plotting_style as ps

def get_timecourse(results_dir, pattern):
    files = os.listdir(results_dir)
    filenames = [f for f in files if
                 pattern in f]
    if len(filenames) < 1:
        raise ValueError("No file found for pattern {}".format(pattern))
    elif len(filenames) > 1:
        raise ValueError(
            "More than one "
            "file found matching {}".format(pattern))
    
    import_file_path = os.path.join(results_dir, filenames[0])

    return pd.read_table(import_file_path, header=0, delim_whitespace=True)    

def plot_wh_timecourse(timecourse,
                       bottleneck,
                       results_dir,
                       wt_col,
                       mut_col,
                       cell_col = None,
                       col_colors = None,
                       detection_threshold = None,
                       detection_limit = None,
                       axis = None,
                       detection_color = None,
                       detection_linestyle = "dotted",
                       E_w = None,
                       t_M = None,
                       transmission_threshold = None,
                       non_trans_alpha = 1,
                       gen_param_file = None,
                       frequencies = True,
                       analytical_frequencies = True,
                       infer_t_M = False):

    if axis is None:
        fig, axis = plt.subplots()

    model_name = os.path.basename(results_dir)

    pattern =  "{}{}_{}".format(model_name, bottleneck, timecourse)

    print("Getting timecourse", pattern)

    data = get_timecourse(
        results_dir,
        pattern)
    
    params = pp.get_params(gen_param_file)
    param_names = {"mu": None,
                   "d_v": None,
                   "R0_wh": None,
                   "r_w": None,
                   "r_m": None,
                   "C_max":None,
                   "k":None,
                   "cross_imm":None,
                   "z_mut": None,
                   "t_M": None,
                   "t_N": None,
                   "emergence_threshold": None}
    
    params = {name: pp.get_param(name,
                                 model_name,
                                 params,
                                 default = default)
              for name, default in param_names.items()}

        
    print("plotting timecourse", pattern)

    if frequencies:
        plot_df_freqs(
            data,
            wt_col,
            mut_col,
            transmission_threshold,
            detection_limit = 1,
            ax = axis,
            col_colors = col_colors,
            non_trans_alpha = non_trans_alpha,
            label = False)
            
        
        if detection_threshold is not None:
            thresholds = np.repeat(detection_threshold, data.time.size)
            axis.axhline(detection_threshold,
                         color=detection_color,
                         linestyle=detection_linestyle)
        if analytical_frequencies:
            transmissible = (data[wt_col] + data[mut_col]) > transmission_threshold
            any_mutant = data[mut_col] > params["emergence_threshold"]
            if np.any(any_mutant):
                row_emerge = data[any_mutant].iloc[0]
                f_0 = (row_emerge["Vm"] /
                       (row_emerge["Vw"] + row_emerge["Vm"]))
                t_emerge = row_emerge["time"]


            else:
                row_emerge, t_emerge = None, None
                
            peak_row = data["Vw"].idxmax()
            t_peak = data["time"].iloc[peak_row]

            trans_f_m = [
                f_m(t_final = t,
                    t_M = t_M,
                    delta = params["k"] * E_w,
                    t_emerge = t_emerge,
                    R0 = params["R0_wh"],
                    d_v = params["d_v"],
                    t_peak = t_peak,
                    f_0 = f_0,
                    mu = params["mu"],
                    ongoing_mutate = True)
                if transmissible[index] else None
                for index, t in enumerate(data.time)]
            non_trans_f_m = [
                f_m(t_final = t,
                    t_M = t_M,
                    delta = params["k"] * E_w,
                    t_emerge = t_emerge,
                    R0 = params["R0_wh"],
                    d_v = params["d_v"],
                    t_peak = t_peak,
                    f_0 = f_0,
                    mu = params["mu"],
                    ongoing_mutate = True)
                if not transmissible[index] and
                any_mutant[index]
                else None
                for index, t in enumerate(data.time)]

            f_m_color = "black"
            analytic_lw = mpl.rcParams['lines.linewidth'] * 0.7
            axis.plot(data.time,
                      non_trans_f_m,
                      color = f_m_color,
                      linestyle = "dotted",
                      alpha = non_trans_alpha * 0.5,
                      lw = analytic_lw)
            axis.plot(data.time,
                      trans_f_m,
                      color = f_m_color,
                      linestyle = "dotted",
                      lw = analytic_lw)

    else:
        plot_df_timecourse(
            data,
            wt_col,
            mut_col,
            transmission_threshold,
            cell_col = cell_col,
            ax = axis,
            col_colors = col_colors,
            non_trans_alpha = non_trans_alpha,
            label = False)

    if t_M is not None:
        axis.axvspan(
            t_M,
            max(np.max(data.time), 50), 
            alpha=0.25 ,
            color="gray")

    axis.grid(b=True,
              which="major")
