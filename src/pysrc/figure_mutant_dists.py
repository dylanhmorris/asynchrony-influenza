#!/usr/bin/env python3

import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import norm
import pandas as pd

import plotting_style as ps
import parse_parameters as pp
from plotting_functions import setup_multipanel, subdivided_hist

from point_of_transmission import poisson_max_sus, neutralize_prob_from_z
from figure_kernel_shift_repl import plot_kernel_shift_repl



def get_data(drift_data_paths = None,
             constant_data_paths = None,
             both_data_paths = None,
             pt_data_paths = None):

    data = {}

    data['drift'] = pd.concat(
        [pd.read_table(path,
                       delim_whitespace = True)
         for path in drift_data_paths])

    data['constant'] = pd.concat(
        [pd.read_table(path,
                       delim_whitespace = True)
         for path in constant_data_paths])
    
    data['both'] = pd.concat(
        [pd.read_table(path,
                       delim_whitespace = True)
         for path in both_data_paths])

    data['pt'] = pd.concat(
        [pd.read_table(path,
                       delim_whitespace = True)
         for path in pt_data_paths])
    
    return data
    

def plot_antigenic_changes(data = None,
                           axis = None,
                           phen_col = "maj_phen",
                           naive_ind = None,
                           phen_bins = None,
                           phen_bounds = [-0.1, 0.1],
                           mut_mean = 0,
                           mut_sd = 0.02,
                           kernel_alpha = ps.kernel_alpha,
                           kernel_color = 'black',
                           kernel_linestyle = 'solid',
                           kernel_lw = None,
                           dash_capstyle = None,
                           **kwargs):

    if axis is None:
        fig, axis = plt.subplots()

    if naive_ind is None:
        naive_ind = data['donor_hist_choice'].max()

    min_phen, max_phen = phen_bounds

    if phen_bins is None:
        phen_bins = np.linspace(min_phen, max_phen, 21)
        print(phen_bins)

    n_bins = len(phen_bins)
    bin_width = phen_bins[1] - phen_bins[0]
    print(bin_width)
    
    culled = data[np.logical_and(
        data[phen_col] > min_phen,
        data[phen_col] < max_phen)]

    culled['delta_phen'] = culled[phen_col].diff()
    changes = culled[np.abs(culled.delta_phen) > 0]

    
    
    # ignore diff at the start of a new chain (crucial!)
    changes = changes[changes.link_no > 0]

    naive = changes['recip_hist_choice'] > naive_ind - 1

    n_chains = data.rng_seed.unique().size
    print("n chains: {}".format(n_chains))

    mean_waiting_time = changes.link_no.mean()
    print("mean waiting time: ", mean_waiting_time)
          
    drift_counts, bins = np.histogram(
        changes[naive].delta_phen,
        bins = phen_bins)
    
    select_counts, bins = np.histogram(
        changes[np.logical_not(naive)].delta_phen,
        bins = phen_bins)
    
    colors = [ps.drift_color, ps.inocs_color]

    drift_counts = drift_counts / n_chains
    select_counts = select_counts / n_chains

    print(np.sum(drift_counts + select_counts))
    
    subdivided_hist([drift_counts, select_counts],
                    bins,
                    colors = colors,
                    normed = False,
                    axis = axis)
    
    xvals = np.linspace(min(phen_bins),
                        max(phen_bins),
                        100)
    
    yvals = norm.pdf(xvals,
                     loc = mut_mean,
                     scale = mut_sd)
    yvals = yvals * bin_width # scale normal to match frequency hist
    axis.plot(xvals, yvals,
              linestyle = kernel_linestyle,
              alpha = kernel_alpha,
              color = kernel_color,
              lw = kernel_lw,
              dash_capstyle = dash_capstyle)


def plot_net_antigenic_changes(data=None,
                               axis=None,
                               phen_col="maj_phen",
                               link_no_col="link_no",
                               ignore_zeros=True,
                               density=False,
                               **kwargs):
    if axis is None:
        fig, axis = plt.subplots()
    chains = data.groupby('rng_seed')
    n_chains = data['rng_seed'].size()
    root = float(chains.head(1)[phen_col][0])
    net_change = (chains.tail(1)[phen_col] - root)
    if ignore_zeros:
        net_change = net_change[np.abs(net_change) > 0]
    net_change.hist(
        color = "lightblue",
        alpha = 0.8,
        ax = axis,
        edgecolor = "black",
        density = density,
        **kwargs)
    y_vals = axis.get_yticks()
    axis.set_yticklabels(
        ['{:0.3f}'.format(x / n_chains) for x in y_vals])


def plot_change_dists(plots,
                      drift_data_paths,
                      constant_data_paths,
                      both_data_paths,
                      pt_data_paths,
                      params_path,
                      original_phen_color = "black",
                      original_phen_lw = 4,
                      original_phen_linestyle = "dashed",
                      **kwargs):

    data = get_data(drift_data_paths,
                    constant_data_paths,
                    both_data_paths,
                    pt_data_paths)

    params = pp.get_params(params_path)

    mut_mean = float(params.get("DEFAULT_CHAIN_MUT_MEAN"))
    mut_sd = float(params.get("DEFAULT_CHAIN_MUT_SD"))
    
    for condition in ['drift', 'constant', 'both', 'pt']:

        df = data[condition]

        # plot antigenic changes
        change_axis = plots['change-dist-' + condition + '-sim']
        change_color = 'gray' if condition == 'drift' else ps.inocs_color
        phen_bounds = [-0.35, 0.35]
        plot_antigenic_changes(data = df,
                               axis = change_axis,
                               mut_mean = mut_mean,
                               mut_sd = mut_sd,
                               phen_bounds = phen_bounds,
                               **kwargs)
        change_axis.set_xlabel("antigenic change")
        if change_axis.is_first_col():
            change_axis.set_ylabel("frequency")
        change_axis.set_xlim(phen_bounds)
        change_axis.grid(b=True)
        change_axis.axvline(0,
                            lw = original_phen_lw,
                            color = original_phen_color,
                            linestyle = original_phen_linestyle)
        
    plots["change-dist-drift-sim"].set_title(
        "naive population\n")
    plots["change-dist-constant-sim"].set_title(
        "immediate recall response\n")
    plots["change-dist-both-sim"].set_title(
        "mucosal antibodies and\n"
        "immediate recall response")
    plots["change-dist-pt-sim"].set_title(
        "mucosal antibodies and\n"
        "realstic (48h) recall response")


def main(drift_data_paths = None,
         constant_data_paths = None,
         pt_data_paths = None,
         both_data_paths = None,
         param_file = None,
         output_path = None):

    ## figure styling / setup
    width = 18.3
    mpl.rcParams['font.size'] = width / 1.5
    mpl.rcParams['lines.linewidth'] = width / 2.5
    kernel_lw = width / 4.5
    mpl.rcParams['ytick.major.pad']= width / 2
    mpl.rcParams['xtick.major.pad']= width / 2
    mpl.rcParams['legend.fontsize']= width * 0.9
    mpl.rcParams['legend.title_fontsize']= width 
    mpl.rcParams['legend.handlelength']= 2

    height = width / 2
    fig = plt.figure(figsize=(width, height))

    nrows = 2
    ncols = 4
    gs = gridspec.GridSpec(nrows, ncols)

    ## multipanel setup
    change_dist_sim_row = 0
    change_dist_ana_row = 1
    drift_col, const_col, both_col, pt_col = 0, 1, 2, 3
    
    plot_positions = [
        {"name": "change-dist-drift-sim",
         "grid_position": np.s_[change_dist_sim_row, drift_col],
         "sharex": None,
         "sharey": None},
        
        {"name": "change-dist-constant-sim",
         "grid_position": np.s_[change_dist_sim_row, const_col],
         "sharex": "change-dist-drift-sim",
         "sharey": "change-dist-drift-sim"},
        
        {"name": "change-dist-both-sim",
         "grid_position": np.s_[change_dist_sim_row, both_col],
         "sharex": "change-dist-drift-sim",
         "sharey": "change-dist-drift-sim"},
        
        {"name": "change-dist-pt-sim",
         "grid_position": np.s_[change_dist_sim_row, pt_col],
         "sharex": "change-dist-drift-sim",
         "sharey": "change-dist-drift-sim"},

        {"name": "change-dist-drift-ana",
         "grid_position": np.s_[change_dist_ana_row, drift_col],
         "sharex": "change-dist-drift-sim",
         "sharey": None},
                        
        {"name": "change-dist-constant-ana",
         "grid_position": np.s_[change_dist_ana_row, const_col],
         "sharex": "change-dist-drift-ana",
         "sharey": "change-dist-drift-ana"},
                                                        
        {"name": "change-dist-both-ana",
         "grid_position": np.s_[change_dist_ana_row, both_col],
         "sharex": "change-dist-drift-ana",
         "sharey": "change-dist-drift-ana"},
                                                    
        {"name": "change-dist-pt-ana",
         "grid_position": np.s_[change_dist_ana_row, pt_col],
         "sharex": "change-dist-drift-ana",
         "sharey": "change-dist-drift-ana"}
    ]

    letter_loc = (-0.1, 1.15)
    plots = setup_multipanel(fig,
                             plot_positions,
                             letter_loc=letter_loc,
                             gridspec=gs)



    ## parametrization
    params = pp.get_params(param_file)
    b = int(params.get("DEFAULT_CHAIN_BOTTLENECK"))
    v = int(float(params.get("DEFAULT_VB_RATIO")) * b)
    f_mut = float(params.get("DEFAULT_F_MUT"))
    z_wt = float(params.get("DEFAULT_Z_WT"))
    k = float(params.get("CONSTANT_CHAIN_K"))
    mu = float(params.get("DEFAULT_MU"))
    d_v = float(params.get("DEFAULT_D_V"))
    R0 = float(params.get("DEFAULT_R0_WH"))
    mut_sd = float(params.get("DEFAULT_CHAIN_MUT_SD"))
    mut_mu = float(params.get("DEFAULT_CHAIN_MUT_MEAN"))
    
    print("parsed parameter file: v = {}, b = {}".format(v, b))
    print("parsed parameter file: f_mut = {}".format(f_mut))
    print("parsed parameter file: z_wt = {}".format(z_wt))
    
    #####################################
    # simulated kernel shifts
    #####################################
    
    print("plotting simulated kernel shifts...")
    plot_change_dists(plots,
                      drift_data_paths,
                      constant_data_paths,
                      both_data_paths,
                      pt_data_paths,
                      param_file,
                      kernel_lw = kernel_lw,
                      kernel_color = ps.kernel_color,
                      kernel_alpha = ps.kernel_alpha,
                      kernel_linestyle = ps.kernel_linestyle,
                      dash_capstyle = 'round')

    #####################################
    # analytical kernel shifts
    #####################################
    print("plotting analytical kernel shifts...")
    ana_phenotypes = np.linspace(-0.35, 0.35, 101)

    ana_gen_phen = -0.8   # least immune history 
    ana_recip_phen = -0.8 # from the sim model

    kern_escape = 1
    kern_model = 'linear'
    
    print("no immunity...")
    plot_kernel_shift_repl(
        axis = plots['change-dist-drift-ana'],
        k = 0,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = b,
        t_M = 0,
        t_transmit = 2,
        vb_ratio = v / b,
        sd_mut = mut_sd,
        phenotypes = ana_phenotypes,
        generator_phenotype =-99,
        z_homotypic = 0.95,
        recipient_phenotype = -99,
        susceptibility_model = kern_model,
        escape = kern_escape,
        dist_color = ps.inocs_color,
        kernel_color = ps.kernel_color,
        kernel_alpha = ps.kernel_alpha,
        kernel_lw = kernel_lw,
        kernel_linestyle = ps.kernel_linestyle,
        dist_alpha = ps.dist_alpha,
        dash_capstyle = 'round',
        mark_original_phenotype = True)

    print("constant recall response...")
    plot_kernel_shift_repl(
        axis = plots['change-dist-constant-ana'],
        k = k,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = b,
        t_M = 0,
        t_transmit = 2,
        vb_ratio = v / b,
        sd_mut = mut_sd,
        phenotypes = ana_phenotypes,
        generator_phenotype = ana_gen_phen,
        z_homotypic = 0.95,
        recipient_phenotype = -99,
        susceptibility_model = kern_model,
        escape = kern_escape,
        dist_color = ps.inocs_color,
        kernel_color = ps.kernel_color,
        kernel_alpha = ps.kernel_alpha,
        kernel_linestyle = ps.kernel_linestyle,
        kernel_lw = kernel_lw,
        dist_alpha = ps.dist_alpha,
        dash_capstyle = 'round',
        mark_original_phenotype = True)

    print("constant recall response with mucosal antibodies...")
    plot_kernel_shift_repl(
        axis = plots['change-dist-both-ana'],
        k = k,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = b,
        t_M = 0,
        t_transmit = 2,
        vb_ratio = v / b,
        sd_mut = mut_sd,
        phenotypes = ana_phenotypes,
        generator_phenotype = ana_gen_phen,
        z_homotypic = 0.95,
        recipient_phenotype = ana_recip_phen,
        susceptibility_model = kern_model,
        escape = kern_escape,
        dist_color = ps.inocs_color,
        kernel_color = ps.kernel_color,
        kernel_alpha = ps.kernel_alpha,
        kernel_linestyle = ps.kernel_linestyle,
        kernel_lw = kernel_lw,
        dist_alpha = ps.dist_alpha,
        dash_capstyle = 'round',
        mark_original_phenotype = True)
    
    print("realistic recall response with mucosal antibodies...")
    plot_kernel_shift_repl(
        axis = plots['change-dist-pt-ana'],
        k = k,
        mu = mu,
        d_v = d_v,
        R0 = R0,
        bottleneck = b,
        t_M = 2,
        t_transmit = 2,
        vb_ratio = v / b,
        sd_mut = mut_sd,
        phenotypes = ana_phenotypes,
        generator_phenotype = ana_gen_phen,
        z_homotypic = 0.95,
        recipient_phenotype = ana_recip_phen,
        susceptibility_model = kern_model,
        escape = kern_escape,
        dist_color = ps.inocs_color,
        kernel_color = ps.kernel_color,
        kernel_alpha = ps.kernel_alpha,
        kernel_lw = kernel_lw,
        kernel_linestyle = ps.kernel_linestyle,
        dist_alpha = ps.dist_alpha,
        dash_capstyle = 'round',
        mark_original_phenotype = True)


    #####################################
    # plot styling
    #####################################

    for plotname, plot in plots.items():
        if 'ana' in plotname:
            plot.set_xlabel('antigenic change')
        if 'drift-sim' in plotname:
            plot.set_ylabel('frequency')
        if 'drift-ana' in plotname:
            plot.set_ylabel('probability density')
        plot.set_ylim(bottom = 0)
        plot.label_outer()
        
    fig.tight_layout()

    # save
    fig.savefig(output_path)

    return 0


if __name__ == "__main__":
    if len(sys.argv) < 1 + 2 + 4 + 4 + 2:
        print("USAGE: {} "
              "<drift data> <const. imm. data> "
              "<const./IgA data> <var./IgA data> "
              "<param file> <output path> "
              "\n\n".format(sys.argv[0]))
    else:
        main(drift_data_paths = sys.argv[1:5],
             constant_data_paths = sys.argv[5:9],
             both_data_paths = sys.argv[9:13],
             pt_data_paths = sys.argv[13:17],
             param_file=sys.argv[17],
             output_path=sys.argv[18])
