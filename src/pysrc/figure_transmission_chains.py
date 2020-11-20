#!/usr/bin/env python3

######################################################
# file: figure_transmission_chains.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: produces publication ready 
# visualization of the transmission chain model
######################################################

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
import plotting_style as ps
import parse_parameters as pp
import matplotlib.gridspec as gridspec
from collections import OrderedDict
import sys


    
if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("\nUSAGE ./figure_antigenic_trajectories.py <drift data> <constant immunity data> <constant w/ point of transmission data> <point of transmission immunity data> <output_path> <parameters>\n\n")
    else:
        main(drift_data_path=sys.argv[1],
             constant_data_path=sys.argv[2],
             both_data_path=sys.argv[3],
             pt_data_path=sys.argv[4],
             output_path=sys.argv[5],
             params_path=sys.argv[6])
