#!/bin/env python3
#
# analyze_selection_comparison.py
#
# helper functions for analyzing simulations
# of the various within-host models
#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os
from distutils.sysconfig import parse_makefile
import re

def import_data(path):
    dat = pd.read_table(path, delim_whitespace=True, comment="#")
    return dat

def get_vb_data(model_output_dir,
                extension=".txt",
                n_encounters_iga=None):
    """
    get data from a set of sims 
    that vary the vb ratio, fixing
    the bottleneck at 1
    """
    if os.path.isdir(model_output_dir):
        files = os.listdir(model_output_dir)
        dir = model_output_dir
    else:
        files = [os.path.basename(model_output_dir)]
        dir = os.path.dirname(model_output_dir)
    df = pd.DataFrame()

    if n_encounters_iga is None:
        pattern = '([0-9]+)' + extension
        add_v = True
    else:
        pattern = '{}'.format(n_encounters_iga) + extension
        add_v = False
    for file in files:
        if re.search(pattern, file):
            filepath = os.path.join(dir, file)
            print("Reading ", filepath)
            temp_dat = import_data(filepath)
            if add_v:
                pat = re.compile(pattern)
                v = int(pat.findall(file)[0])
                temp_dat['n_encounters_iga'] = v
            df = pd.concat([df, temp_dat], axis=0)
    return df


def get_bottleneck_data(model_output_dir,
                         extension=".txt"):
    """
    get data from a set of sims that
    vary the bottleneck 
    """
    if os.path.isdir(model_output_dir):
        files = os.listdir(model_output_dir)
        dir = model_output_dir
    else:
        files = [os.path.basename(model_output_dir)]
        dir = os.path.dirname(model_output_dir)
    df = pd.DataFrame()

    pattern = '([0-9]+)' + extension
    for file in files:
        if re.search(pattern, file):
            filepath = os.path.join(dir, file)
            print("Reading ", filepath)
            temp_dat = import_data(filepath)
            pat = re.compile(pattern)
            bottleneck = int(pat.findall(file)[0])
            temp_dat['bottleneck'] = bottleneck
            df = pd.concat([df, temp_dat], axis=0)
    return df


def get_sensitivity_paramsets(model_output_dir,
                  extension='.txt',
                  bottleneck=None,
                  paramset_string='_paramsets'):
    """
    get the record of random parameter 
    sets generated during the sensitivity
    analysis (which deterministically varies
    the bottleneck)
    """
    if os.path.isdir(model_output_dir):
        files = os.listdir(model_output_dir)
        dir = model_output_dir
    else:
        files = [os.path.basename(model_output_dir)]
        dir = os.path.dirname(model_output_dir)
    df = pd.DataFrame()

    if bottleneck is None:
        pattern = '([0-9]+)' + paramset_string + extension
        add_bn = True
    else:
        pattern = '{}'.format(bottleneck) + paramset_string  + extension
        add_bn = False
        
    for file in files:
        if re.search(pattern, file):
            filepath = os.path.join(dir, file)
            print("Reading ", filepath)
            temp_dat = import_data(filepath)
            if add_bn:
                pat = re.compile(pattern)
                bottleneck = int(pat.findall(file)[0])
            temp_dat['bottleneck'] = bottleneck
            df = pd.concat([df, temp_dat], axis=0)
    return df


def get_mean_fmut(df):
    return df["p_mut_inoc"].groupby(
        df["bottleneck"]).mean()


def get_p_select_per_infection(df, freq):
    subset = df[df.transmissible]
    if subset.size() > 0:
        print(subset.head())
        denom = subset.groupby("bottleneck").size()
        selection_prob = (subset.groupby(["bottleneck"]
        ).mean_mutant_freq.apply(
            lambda x: (x > freq).sum())
                          / denom)
    else:
        print("Warning: no transmissible cases! Maybe run more simulations")
        selection_prob = 0
    return selection_prob


def get_selection_probs(df, freq):
    if df.size > 0:

        selection_prob = (df.groupby("bottleneck").mean_mutant_freq.apply(
            lambda x: (x > freq).sum())
                          / df.groupby("bottleneck").size())
    
        repl_prob = (
            df[np.logical_not(df.is_inoc)].groupby("bottleneck"
            ).mean_mutant_freq.apply(lambda x: (x > freq).sum())
            / df.groupby("bottleneck").size())
    
        inocs_prob = selection_prob - repl_prob
    else:
        print("Warning: no cases! Maybe run more simulations")
        return None
    return (selection_prob, repl_prob, inocs_prob)
