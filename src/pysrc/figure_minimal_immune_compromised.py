#!/usr/bin/env python3

##################################################
# filename: figure_minimal_immune_compromised.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: generates figure showing that
# replication selection occurs in immune compromised
# hosts
##################################################

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

import plotting_style as ps
from plotting_functions import setup_multipanel
from wh_popgen import f_m

def setup_figure():
    width = 18.3
    height = (1/3) * width

    lineweight = width / 3
    mpl.rcParams['lines.linewidth'] = lineweight
    mpl.rcParams['axes.formatter.limits'] = (-3, 6)

    # set up multi-panel figure
    fig = plt.figure(figsize=(width, height))
    gs = mpl.gridspec.GridSpec(1, 3)


    short_col, long_col, exp_long_col = 0, 1, 2
    
    plot_positions = [
        {"name": "typical",
         "grid_position": np.s_[0, short_col],
         "sharex": None,
         "sharey": None},
        
        {"name": "prolonged",
         "grid_position": np.s_[0, long_col],
         "sharex": None,
         "sharey":  "typical"},
                
        {"name": "prolonged-experienced",
         "grid_position": np.s_[0, exp_long_col],
         "sharex": "prolonged",
         "sharey": "typical"}
    ]
    
    plots = setup_multipanel(fig,
                             plot_positions,
                             gridspec=gs)

    return (fig, plots)


def main(bottleneck,
         gen_param_file,
         outfile):
    figsize = (10, 5)
    legend_handlelength = 2
    legend_loc = "best"
    tick_step = 5
    detection_linestyle="dashed"

    typical_title = "Typical, short-lived"
    prolonged_title = "Prolonged"
    prolonged_experienced_title = "Prolonged, experienced"

    fig, plots = setup_figure()
    
    print("Plotting infections...")
    handles = []
    print("Plotting competent...")

    times = np.linspace(0, 100, 10000)

    t_M_naive = 6
    t_M_experienced = 2
    delta_typical = 6
    delta_prolonged = 0.25

    R0 = 5
    d_v = 4
    mu = .33e-5
    
    t_emerge = 0.75
    f_0 = mu
    
    f_m_typical = np.array([
        f_m(t_final = time,
            t_M = t_M_naive,
            t_emerge = t_emerge,
            delta = delta_typical,
            R0 = R0,
            d_v = d_v,
            f_0 = f_0,
            mu = mu,
            t_peak = 2,
            ongoing_mutate = False)
        for time  in times])
    
    f_w_typical = 1 - f_m_typical

    axis = plots['typical'].plot(
        times, f_w_typical,
        color = ps.wt_color,
        label = "old variant")
    axis = plots['typical'].plot(
        times, f_m_typical,
        color = ps.mut_color,
        label = "new variant")
        
    print("Plotting prolonged...")
    f_m_prolonged = np.array([
        f_m(t_final = time,
            t_M = t_M_naive,
            t_emerge = t_emerge,
            delta = delta_prolonged,
            R0 = R0,
            d_v = d_v,
            f_0 = f_0,
            mu = mu,
            t_peak = 2,
            ongoing_mutate = False)
        for time  in times])
    
    f_w_prolonged = 1 - f_m_prolonged

    
    axis = plots['prolonged'].plot(
        times, f_w_prolonged,
        color = ps.wt_color)
    axis = plots['prolonged'].plot(
        times, f_m_prolonged,
        color = ps.mut_color)

    print("Plotting prolonged/experienced...")
    f_m_experienced = np.array([
        f_m(t_final = time,
            t_M = t_M_experienced,
            t_emerge = t_emerge,
            delta = delta_prolonged,
            R0 = R0,
            d_v = d_v,
            f_0 = f_0,
            mu = mu,
            t_peak = 2,
            ongoing_mutate = False)
        for time  in times])
    
    f_w_experienced = 1 - f_m_experienced
    
    axis = plots['prolonged-experienced'].plot(
        times, f_w_experienced,
        color = ps.wt_color)
    axis = plots['prolonged-experienced'].plot(
        times, f_m_experienced,
        color = ps.mut_color)

    print("Styling...")
    plots['prolonged'].set_xlabel(r'Time (days)')
    plots['typical'].set_xlim([-0.5, 10])
    plots['prolonged'].set_xlim([0, 100])
    plots['typical'].set_ylabel('variant frequency')
    plots['typical'].set_title(typical_title)
    plots['prolonged'].set_title(prolonged_title)
    plots['prolonged-experienced'].set_title(prolonged_experienced_title)

    plots['typical'].set_xticks(np.arange(0, 12, tick_step))
    
    for axis in plots.values():
        axis.label_outer()
        axis.tick_params(axis='both', which='major')
        axis.tick_params(axis='both', which='minor')
        x_start, x_end = axis.get_xlim()
        axis.xaxis.set_ticks(np.arange(0, x_end, tick_step),
                             minor=True)
        axis.grid(which="minor")
        axis.grid(b=True)
        axis.set_xlim([0, x_end])
        axis.set_yscale("log")
    
    leg = plots['typical'].legend()
    
    for legobj in leg.legendHandles:
        legobj.set_linewidth(ps.standard_lineweight)

    fig.tight_layout()
    fig.savefig(outfile)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("\nUSAGE: ./figure_minimal_immune_compromised.py "
              "bottleneck "
              "gen_param_file "
              "outfile\n")
    else:
        main(sys.argv[1],
             sys.argv[2],
             sys.argv[3])
