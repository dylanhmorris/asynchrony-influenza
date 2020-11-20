#!/usr/bin/env python3

# Plots the analytic escape probability

from analytic_escape_prob import get_basic_prob, get_q
from multi_class_final_size import multi_class_final_size
import plotting_style as ps
import sys
import numpy as np
import matplotlib.pyplot as plt
from figure_db_tools import output_figure_record
import os

def main(outpath, N, h):
    axis_fontsize = ps.standard_axis_fontsize
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_analytical(N, h, axis=ax)
    caption_templ = (
        "Probability of at least one mutant selective "
        "event during an epidemic as a function of the "
        "initial fraction immune in the population. "
        "Mutant assumed to offer complete escape. Parameters: "
        "$h$ = {}, $N$ = {}")

    ax.set_title("epidemic-level selection")

    fig_caption = caption_templ.format(h, N)
    fig_name = os.path.splitext(os.path.basename(outpath))[0]
    fig.savefig(outpath)

def plot_q(axis=None):
    if not axis:
        fig, axis = plt.subplots()

    axis_fontsize = ps.standard_axis_fontsize
    axis_ticksize = ps.standard_axis_ticksize

    R0s = [1.1, 1.4, 1.7, 2]
    s0s = np.arange(0, 1.01, 0.01)[::-1]
    mylines = []
    for R0 in R0s:
        xs = 1 - s0s
        approx_probs = [get_q(s0, R0) for s0 in s0s]
        mylines += axis.plot(xs, approx_probs, label=str(R0))

def plot_analytical(N, h, axis=None,
                    cmap=plt.cm.Reds,
                    per_capita=True):
    if not axis:
        fig, axis = plt.subplots()

    mylines = []

    s0s = np.arange(0, 1.01, 0.01)[::-1]
    R0s = [1.2, 1.5, 1.8]
    line_ids = np.linspace(0.15, 0.7, len(R0s))

    lineweight = ps.standard_lineweight
    axis_fontsize = ps.standard_axis_fontsize
    axis_ticksize = ps.standard_axis_ticksize

    for R0, line_id in zip(R0s, line_ids):
        q_naive = get_q(1, R0)
        xs = 1 - s0s
        approx_probs = [get_basic_prob(s0, R0, h, N)/N for s0 in s0s]
        mylines += axis.plot(xs, approx_probs, label=str(R0),
                             color=cmap(line_id))
        axis.axvline(x=q_naive,
                     color=cmap(line_id),
                     linestyle="dashed")
    axis.set_xlabel("level of immunity")
    if per_capita:
        axis.set_ylabel("probability per capita")
    else:
        axis.set_ylabel("probability")
    axis.set_xlim([0, 1])
    axis.set_ylim(ymin=0)
    legend = axis.legend(handles=mylines,
                         frameon=True,
                         fancybox=True,
                         framealpha=1,
                         labelspacing=0.25,
                         handlelength=0.5)
    legend.set_title("$\mathcal{R}_0$")
    axis.grid()


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Please supply output path, N, and h")
    else:
        main(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]))
