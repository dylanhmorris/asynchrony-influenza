#!/usr/bin/env python3

##############################################
# name: plotting_style.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: fix general plotting style
# for the project
##############################################


import matplotlib as mpl
import matplotlib.pyplot as plt

plt.style.use("seaborn-white")

mpl.rcParams['axes.labelsize'] = "large"
mpl.rcParams['xtick.labelsize'] = "large"
mpl.rcParams['ytick.labelsize'] = "large"
mpl.rcParams['axes.titlesize'] = "x-large"
mpl.rcParams['axes.formatter.use_mathtext'] = True
mpl.rcParams['axes.formatter.limits'] = ((-3, 3))
mpl.rcParams['axes.grid'] = True
mpl.rcParams['legend.frameon'] = True
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.framealpha'] = 1
mpl.rcParams['legend.fontsize'] = "medium"
mpl.rcParams['lines.linewidth'] = 5


mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{amssymb}']

mpl.rcParams['text.usetex'] = True
mpl.rcParams['axes.titleweight']="bold"
mpl.rcParams['mathtext.fontset'] = 'stix'

journ_width = 18.3
mpl.rcParams['ytick.major.pad'] =  journ_width / 2
mpl.rcParams['xtick.major.pad'] = journ_width / 2

mpl.rcParams['font.size'] = journ_width

###########################
## color scheme
###########################

wt_color = "#4fbfd1"
mut_color = "#db5851"
cell_color = "#7aa081"
no_inf_color = "#f7f4f2"
mixed_inoc_color = "#600087"

antigenic_site_color = "#6c97ab"
antigenic_ridge_color = "#ff0000"
antigenic_ridge_color = "#ff0000"
ridge_vax_color = "black"
generic_substitution_color = "#3a464d"
no_substitution_color = "#f8fae1"

f_color = "purple"
a_color = "brown"
trans_threshold_color = "black"
detection_color = "black"
immune_active_color = "gray"

inocs_color = "#14b200"
repl_color = "#9c40ed"
both_color = "#c2000a"
drift_color= "#848587"

old_phen_color = "#729CDF"
new_phen_color = "#FDD34A"
null_phen_color = "#8A9197"

## kernel and distribution
## styling for distribution shift
## figures
kernel_color = "black"
kernel_alpha = 0.6
dist_alpha = 1
kernel_linestyle = 'dotted'

########################
## general params
########################
standard_lineweight = 5
standard_axis_fontsize = 30
standard_axis_ticksize = 25

letter_loc = (-0.05, 1.15)
letter_size = "xx-large"


#######################################
# how to render parameters as LaTeX
#######################################
parm_display_names = {
    "R0_WH": "$\mathcal{R}_{0}$",
    "R_W": "$r_w$",
    "MU": "$\mu$",
    "D_V": "$d_V$",
    "K": "$k$",
    "CROSS_IMM": "$c$",
    "T_M": "$t_M$",
    "TD50": "$V_{50}$",
    "T_N": "$t_N$",
    "C_MAX": "$C_{max}$",
    "Z_WT": "$z_{w}$",
    "MUT_WT_NEUT_RATIO": "$\kappa_m / \kappa_w$",
    "Z_MUT_Z_WT_RATIO": "$z_{m} / z_{w}$",
    "TRANS_THRESHOLD": "trans. threshold",
    "VB_RATIO": "$v/b$",
    "bottleneck": "bottleneck"}
