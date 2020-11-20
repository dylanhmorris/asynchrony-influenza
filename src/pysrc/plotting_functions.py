import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_cols
import matplotlib.image as mpimg
import plotting_style
import pandas as pd
from matplotlib.ticker import ScalarFormatter, LogFormatterMathtext
from textwrap import wrap
import plotting_style as ps

def add_bounding_subplot(figure, position=None):
    if position is None:
        position = 111
    ax = figure.add_subplot(position)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.grid(b=False)
    ax.patch.set_alpha(0)
    ax.tick_params(labelcolor='w',
                   grid_alpha=0,
                   top=False,
                   bottom=False,
                   left=False,
                   right=False)
    for label in ax.get_xticklabels():
        label = ''
    for label in ax.get_yticklabels():
        label = ''
    ax.set_zorder(0)
    return ax


def latex_pow_10(float):
    exponent = np.log10(float)
    the_exp = int(exponent)
    if abs(the_exp - exponent) > 0.1:
        raise ValueError("This number does not appear to be close to "
                         "a power of 10")
    else:
        return r"$10^{{{0}}}$".format(the_exp)


def scale_color_brightness(color, factor, darken=True):
    rgb = mpl_cols.to_rgb(color)
    if factor > 1 or factor < 0:
        raise ValueError("Invalid scaling factor; "
                         "must be between 0 and 1")
    if darken:
        shade_factor = factor
        new_rgb = [current * (1 - shade_factor)
                   for current in rgb]
    else:
        tint_factor = factor
        new_rgb = [current + (255 - current) * tint_factor
                   for current in rgb]
    return new_rgb


def plot_df_timecourse(dataframe,
                       wt_col,
                       mut_col,
                       transmission_threshold,
                       col_colors = None,
                       time_col = "time",
                       cell_col=None,
                       ax = None,
                       label=True,
                       display_names = None,
                       non_trans_alpha = 1,
                       **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    if display_names is None:
        display_names = ["wild-type", "mutant", "target cells"]

    if not col_colors:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        col_colors = [colors[ind] for ind in range(3)]

    col_names = [wt_col, mut_col]

    total_pop = dataframe[wt_col] + dataframe[mut_col]
    transmissible = total_pop > transmission_threshold
    
    handles = []

    if cell_col is not None:
        ax.plot(dataframe[time_col],
                dataframe[cell_col],
                color=col_colors[2],
                label=display_names[2],
                **kwargs)

    for column, display_name, color in zip(col_names,
                                           display_names,
                                           col_colors):
        trans_vals = [item if transmissible[index] else None for
                       index, item in enumerate(dataframe[column])]
        non_trans_vals = [item if not transmissible[index] else None for
                           index, item in enumerate(dataframe[column])]

        handles += ax.plot(dataframe[time_col],
                           trans_vals,
                           label=display_name,
                           color=color,
                           **kwargs)
        ax.plot(dataframe[time_col],
                non_trans_vals,
                color=color,
                alpha=non_trans_alpha,
                **kwargs)
        ax.set_yscale('log')

        
    if label:
        plt.xlabel(r'time (days)')
        plt.ylabel(r'virions')
        plt.legend(handles=handles)
    return handles


def plot_df_freqs(dataframe,
                  wt_col,
                  mut_col,
                  transmission_threshold,
                  display_names = None,
                  col_colors = None,
                  time_col = "time",
                  detection_limit = 1,
                  ax = None,
                  label = True,
                  non_trans_alpha = 1,
                  **kwargs):
    if not ax:
        fig, ax = plt.subplots()

    if display_names is None:
        display_names = ["wild-type", "mutant"]

    if col_colors is None:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        col_colors = [colors[ind] for ind in range(3)]

    col_names = [wt_col, mut_col]

    total_pop = dataframe[wt_col] + dataframe[mut_col]
    transmissible = total_pop > transmission_threshold
    detectable = total_pop > detection_limit

    handles = []
    for column, display_name, color in zip(col_names,
                                           display_names,
                                           col_colors):
        frequencies = dataframe[column]/dataframe[col_names].sum(axis=1)
        trans_freqs = [item if detectable[index] and
                       transmissible[index] else None for
                       index, item in enumerate(frequencies)]
        non_trans_freqs = [item if not transmissible[index] and
                           detectable[index] else None for
                           index, item in enumerate(frequencies)]
        
        ax.plot(dataframe[time_col],
                non_trans_freqs,
                color=color,
                alpha=non_trans_alpha,
                **kwargs)
        handles += ax.plot(dataframe[time_col],
                           trans_freqs,
                           label=display_name,
                           color=color,
                           **kwargs)
        ax.set_yscale('log')
        
    if label:
        plt.xlabel(r'time (days)')
        plt.ylabel(r'frequency')
        plt.legend(handles=handles)
    return handles



def plot_df_detectable_level(df,
                             wt_col,
                             mut_col,
                             time_col = "time",
                             detection_threshold = 0.01,
                             detectable_population = 1e5,
                             detection_label="detection threshold",
                             linestyle="solid",
                             linewidth = 6,
                             linecolor = "black",
                             alpha = 1,
                             ax = None):
    if not ax:
        fig, ax = plt.subplots()

    handles = []
    times = df[time_col]
    total_pop = df[mut_col] + df[wt_col]

    detectable_level = total_pop * detection_threshold
    threshold = total_pop > detectable_population
    
    handles += ax.plot(times[threshold],
                       detectable_level[threshold],
                       label=detection_label,
                       linestyle=linestyle,
                       linewidth=linewidth,
                       color=linecolor,
                       alpha=alpha)
    return handles



def plot_infection_course(system, wt_col, mut_col, label=True,
                          wt_label="wild-type", mut_label="mutant",
                          mut_color=None, wt_color=None,
                          ax=None, **kwargs):
    if not ax:
        fig, ax = plt.subplots()
    handles = []
    wt_ts = [item[wt_col] for item in system.array_ts]
    mut_ts = [item[mut_col] for item in system.array_ts]
    ax.set_yscale('log')

    handles += ax.plot(system.time_ts,
                       wt_ts,
                       label=wt_label,
                       color=wt_color,
                       **kwargs)
    handles += ax.plot(system.time_ts,
                       mut_ts,
                       label=mut_label,
                       color=mut_color,
                       **kwargs)
    
    if label:
        plt.xlabel(r'Time (days)')
        plt.ylabel(r'Free virions')
        plt.legend(handles=handles)
    
    return handles


def plot_detectable_level(system, wt_col, mut_col,
                          detection_threshold=0.01,
                          detectable_population = 1e5,
                          detection_label="detection threshold",
                          linestyle="solid",
                          linewidth=6,
                          linecolor="black",
                          alpha = 1,
                          ax=None):
    if not ax:
        fig, ax = plt.subplots()
    ax.set_yscale('log')

    handles = []
    times = np.array(system.time_ts)
    total_pop = np.array([item[mut_col] + item[wt_col]
                           for item in system.array_ts])

    detectable_level = total_pop * detection_threshold
    threshold = total_pop > detectable_population
    
    handles += ax.plot(times[threshold],
                       detectable_level[threshold],
                       label=detection_label,
                       linestyle=linestyle,
                       linewidth=linewidth,
                       color=linecolor,
                       alpha=alpha)
    return handles


def stacked_bar(groups,
                trials,
                successes,
                xlab="",
                ylab="Frequency",
                success_label="Mutant seen",
                failure_label="Mutant not seen",
                title=None,
                legend=True,
                failure_color=None,
                success_color=None,
                failure_alpha=0.05,
                figsize=(35, 8),
                normed=False,
                ax=None,
                **kwargs):
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)
    width = 0.8
    if normed:
        successes = successes / trials
        trials = [1 for _ in successes]
    ind = [x for x, _ in enumerate(groups)]
    ax.bar(ind, trials,
           width=width,
           color=failure_color,
           label=failure_label,
           alpha=failure_alpha,
           **kwargs)

    ax.bar(ind, successes, width=width, color=success_color, label=success_label, **kwargs)

    ax.set_xticks(ind)
    ax.set_xticklabels(groups, fontsize = 10)
    ax.set_ylabel(ylab, fontsize=15)
    ax.set_xlabel(xlab, fontsize=15)
    if legend:
        plt.legend(fontsize=14)
    if title:
        ax.set_title(title, fontsize=16)


def get_position(position_dict,
                 n_cols):
    pos = (position_dict["row"] - 1) * n_cols + position_dict["col"]
    return pos


def setup_multipanel(fig, plot_positions,
                     gridspec = None,
                     add_letters = True,
                     letter_loc = None,
                     verbose = False,
                     letters = None,
                     upper = True):

    if gridspec is None:
        n_rows = np.max([plot_dict["row"] for plot_dict in plot_positions])
        n_cols = np.max([plot_dict["col"] for plot_dict in plot_positions])
        
    plots = {}

    if letters is None:
        letters = ["a", "b", "c", "d", "e", "f",
                   "g", "h", "i", "j", "k", "l",
                   "m", "n", "o", "p", "q", "r",
                   "s", "t", "u", "v", "w", "x",
                   "y", "z", "aa", "bb", "cc"]

    if upper:
        letters = [let.upper() for let in letters]

    for plot_no, plot_dict in enumerate(plot_positions):

        sharex = plot_dict["sharex"]
        sharey = plot_dict["sharey"]

        if sharex is not None:
            plot_sharex = plots[sharex]
        else:
            plot_sharex = sharex
            
        if sharey is not None:
            plot_sharey = plots[sharey]
        else:
            plot_sharey = sharey

        if gridspec is None:
            position_id = get_position(plot_dict,
                                       n_cols)
            if verbose:
                print("Adding plot {} at position {}"
                      "".format(plot_dict['name'],
                                    position_id))
            ax = fig.add_subplot(n_rows,
                                 n_cols,
                                 position_id,
                                 sharex=plot_sharex,
                                 sharey=plot_sharey)
        else:
            grid_pos = plot_dict["grid_position"]
            
            if verbose:
                print("Adding plot {} at position {}"
                      "".format(plot_dict['name'],
                                    grid_pos))

            ax = fig.add_subplot(gridspec[grid_pos],
                                 sharex=plot_sharex,
                                 sharey=plot_sharey)
        if add_letters:
            letter_locations = plot_dict.get("letter_loc", letter_loc)
            if letter_locations is None:
                letter_x, letter_y = ps.letter_loc
                # handle spanny letters equivalently
                # by default
                if gridspec is not None:
                    (_, _,
                     row_start,
                     row_stop,
                     col_start,
                     col_stop) = ax.get_subplotspec().get_rows_columns()
                    colspan = 1 + col_stop - col_start
                    rowspan = 1 + row_start - row_stop
                    x_anchor = (
                        0 if abs(1 - letter_x) > abs(0 - letter_x) else 1)
                    y_anchor = (
                        0 if abs(1 - letter_y) > abs(0 - letter_y) else 1)
                    l_x_rel = letter_x - x_anchor
                    l_y_rel = letter_y - y_anchor
                    letter_x = x_anchor + l_x_rel / colspan
                    letter_y = y_anchor + l_y_rel / rowspan                  
            elif type(letter_locations) is tuple:
                letter_x, letter_y = letter_locations
            if plot_dict.get("include_letter", True):
                ax.text(letter_x, letter_y, letters[plot_no],
                        transform=ax.transAxes,
                        fontsize=ps.letter_size,
                        fontweight='bold', va='top')
        
        plots[plot_dict["name"]] = ax

    return plots


def subdivided_hist(vector_of_counts,
                    bins = None,
                    normed = True,
                    axis = None,
                    edgecolor = 'k',
                    cmap = None,
                    colors = None,
                    labels = None,
                    legend = False,
                    **kwargs):
    """
    Plots a histogram with bars 
    subdivided to show within-bin
    variation according to some 
    other provided variable
    """

    if axis is None:
        fig, axis = plt.subplots()
    
    if bins is None:
        bins = np.arange(len(vector_of_counts[0]) + 1)
    n_bins = len(bins)

    if labels is None:
        labels = [None] * len(vector_of_counts)
    
    if not all([len(count_vec) == n_bins - 1
                for count_vec in vector_of_counts]):
        raise ValueError("all vectors of counts must "
                         "have the same number of bins, and it must "
                         "match the provided bins vector "
                         "(which also includes the right endpoint), "
                         "if one is provided\n")
    if n_bins < 2:
        raise ValueError("histogram must have at least two bins")
    
    width = (bins[1] - bins[0])
    inds = np.array([item + width / 2 for item in bins][:-1])
    
    total_count = np.sum(vector_of_counts)

    if normed:
        to_plot = [np.array(vec) / total_count
                   for vec in vector_of_counts]
    else:
        to_plot = [np.array(vec) for vec in
                   vector_of_counts]
    
    for ind, count_vec in enumerate(to_plot):
        if cmap is not None:
            col = cmap(ind)
        elif colors is not None:
            col = colors[ind]
        else:
            col = None

        bottom = None
        if ind > 0:
            bottom = np.sum(to_plot[0:ind],
                            axis=0)
            
        axis.bar(inds,
                 count_vec,
                 width,
                 color = col,
                 bottom = bottom,
                 edgecolor = edgecolor,
                 label = labels[ind],
                 **kwargs)
    if legend:
        axis.legend()

def plot_image_from_path(path,
                         axis = None,
                         retain_axis = False):
    """
    Convenience wrapper 
    for reading in and
    plotting images
    from an image path
    """
    if axis is None:
        fig, axis = plt.subplots()

    img = mpimg.imread(path)
    axis.imshow(img)

    if not retain_axis:
        axis.axis('off')
