#!/usr/bin/env python3

##############################################
# name: process_wh_data.py
# author: Dylan Morris <dhmorris@princeton.edu>
# description: clean data from other studies
# for reanalysis, including performing some
# data integrity checks
##
## studies included:
##
## 1) Debbink, K et al. Vaccination has minimal impact on the
## intrahost diversity of H3N2 influenza viruses.
## PLoS Pathogens. 2017.
## DOI:https://doi.org/10.1371/journal.ppat.1006194
##
## 2) McCrone, JT et al. Stochastic processes constrain
## the within and between host evolution of influenza virus
## eLife. 2018. 
## DOI: https://doi.org/10.7554/eLife.35962.001
##
## All data used is public and
## linked from the published articles above
##############################################

## imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os


## needed constants and values 

antigenic_sites = {
    "A": [122, 124, 126, 130, 131,
          132, 133, 135, 137, 138,
          140, 142, 143, 144, 145,
          146, 150, 152, 168],
    
    "B": [128, 129, 155, 156, 157,
          158, 159, 160, 163, 164,
          165, 186, 187, 188, 189,
          190, 192, 193, 194, 196,
          197, 198],

    "C": [44, 45, 46, 47, 48, 50,
          51, 53, 54, 273, 275,
          276, 278, 279, 280, 294,
          297, 299, 300, 304, 305,
          307, 308, 309, 310, 311,
          312], 

    "D": [96, 102, 103, 117, 121,
          167, 170, 171, 172, 173,
          174, 175, 176, 177, 179,
          182, 201, 203, 207, 208,
          209, 212, 213, 214, 215,
          216, 217, 218, 219, 226,
          227, 228, 229, 230, 238,
          240, 242, 246, 247, 248], 

    "E": [57, 59, 62, 63, 67, 75,
          78, 80, 81, 82, 83, 86,
          87, 88, 91, 92, 94, 109,
          260, 261, 262, 265]
}

all_antigenic_sites = pd.Series([item for sublist in
                                 antigenic_sites.values()
                                 for item in sublist])
                      
antigenic_ridge = pd.Series(
    [145, 155, 156, 158,
     159, 189, 193])


mc_dtypes = {
    'Id': 'str',
    'LAURING_ID': 'str',
    'dup': 'str',
    'SPECID': 'str',
    'ENROLLID': 'str',
    'ENROLLID.1': 'str'}

column_spec = ['Id',
               'season',
               'vaccination_status',
               'vaccination_matched',
               'ref',
               'var',
               'aa_pos',
               'aa_consensus',
               'aa_var',
               'H3_pos',
               'is_aa_substitution',
               'freq_var',
               'antigenic_site',
               'is_antigenic_site',
               'is_antigenic_ridge']
    



#################################################
## useful functions for applying to dataframes
#################################################

def is_site(position_array):
    """
    Given a position in H3N2 HA,
    assess whether it lies in a classically-
    defined antigenic site per Wiley, DC et 
    al. 1981, "Structural identification of the 
    antibody-binding sites of Hong Kong influenza 
    haemagglutinin and their involvement 
    in antigenic variation", Nature.
    DOI: https://doi.org/10.1038/289373a0
    """
    return position_array.isin(all_antigenic_sites)

def is_ridge(position_array):
    """
    Given a position in H3N2 HA,
    assess whether it lies in the 
    "antigenic ridge" identified
    by Koel et al. 2013 "Substitutions Near 
    the Receptor Binding Site Determine 
    Major Antigenic Change During 
    Influenza Virus Evolution. Science.
    DOI: https://doi.org/10.1126/science.1244730
    """

    return position_array.isin(antigenic_ridge)


def which_site(position):
    """
    Given a position in H3N2 HA,
    return the antigenic site (Wiley
    et al. 1981) in which it lies,
    if any.
    """
    result = ""
    if not np.isnan(position):
        for site, positions in antigenic_sites.items():
            if position in positions:
                result = site
    return result

def get_sites(position_array):
    """
    vectorized version of which_site()
    """
    return position_array.apply(which_site)

def compute_mutants(dat):
    """
    Compute whether a variant 
    is a mutant relative to the reference
    genome; this may be the major allele.
    """
    dat["mut_major"] = dat["ref"] == dat["var"]
    dat["freq_mut"] = dat.apply(
        lambda x: 1 - x["freq.var"] if x["mut_major"]
        else x["freq.var"],
        axis = 1)
    dat["phen_mut"] = dat.apply(
        lambda x: x["aa_consensus"] if
        x["mut_major"] else x["aa_var"],
        axis = 1)
    dat["phen_ref"] = dat.apply(
        lambda x: x["aa_var"] if
        x["mut_major"] else x["aa_consensus"],
        axis = 1)
    return dat

def compute_antigenic(dat, pos_column = "H3_pos"):
    """
    Convenience function to do antigenic analysis 
    of a dataframe given a column with H3N2 HA positions
    """
    dat["antigenic_site"] = get_sites(dat[pos_column])
    dat["is_antigenic_site"] = is_site(dat[pos_column])
    dat["is_antigenic_ridge"] = is_ridge(dat[pos_column])
    return dat

def fill_infection_nas(dat):
    """
    Convenience function to replace NAs after 
    joining a dataframe with nonsynonymous 
    substitutions to a larger one representing 
    all sequenced infections, so that the infections 
    without nonsynonymous substitutions are correctly
    enumerated and labeled
    """
    replacements = [
        [0.0, ["freq_var",
               "freq_mut"]],
        
        ["", ["phen_mut",
              "phen_ref",
              "antigenic_site"]],

        [False, ["is_aa_substitution",
                 "is_antigenic_site",
                 "is_antigenic_ridge"]]
    ]

    for replacement, cols in replacements:
        for col in cols:
            if col in dat.columns:
                dat[col] = dat[col].fillna(replacement)
            pass
        pass
    return dat


def vaccination_class(patient):
    """
    convenient function that 
    returns a string vaccination class column 
    a dataframe with a vaccination_status 
    boolean column and a vaccination_match 
    boolean column
    """
    if (patient["vaccination_status"] and
        patient["vaccination_matched"]):
        result = "matched"
    elif (patient["vaccination_status"] and not
          patient["vaccination_matched"]):
        result = "not matched"
    else:
        result = "unvaccinated"
    return result 

def check_data_integrity(dat):
    """
    performs integrity checks for merged and 
    cleaned data
    """

    ## number with mutant frequency > 0
    ## should be same as number with
    ## AA substitution
    freq_check = (dat["is_aa_substitution"].sum() ==
                  (dat["freq_var"] > 0).sum())
    
    if not freq_check:
        raise ValueError("Failed data integrity check: "
                         "mutant frequencies > 0 do not "
                         "match number of rows flagged "
                         "as AA substitutions\n")

    ## all antigenic site substitutions must
    ## also be AA substitutions
    all_site_are_subst = (
        not any(dat.is_antigenic_site &
                (~dat.is_aa_substitution)))

    if not all_site_are_subst:
        raise ValueError("Failed data integrity check: "
                         "at least one row flagged as an "
                         "antigenic site substitution is "
                         "not flagged as an AA substitution\n")

    ## all antigenic ridge substitutions must
    ## also be antigenic site substitutions
    all_ridge_are_site = (
        not any(dat.is_antigenic_ridge &
                (~dat.is_antigenic_site)))

    if not all_ridge_are_site:
        raise ValueError("Failed data integrity check: "
                         "at least one row flagged as an "
                         "antigenic ridge substitution is "
                         "not flagged as an antigenic "
                         "site substitution\n")

    return dat
    
def get_debbink(debbink_data_dir = None,
                suffix = None):
    """
    Read in data from Debbink et al 2017
    and concatenate multiple years in a 
    way that preserves integrity
    """

    debbink_years = [
        "2004-2005",
        "2005-2006",
        "2007-2008"]
        
    debbink_files = [
        "{}.{}.csv".format(year, suffix)
        for year in debbink_years]
    all_dat = pd.DataFrame()
    for season, file in zip(debbink_years, debbink_files):
        path = os.path.join(debbink_data_dir,
                            file)
        dat = pd.read_csv(path)
        dat["season"] = season
        all_dat = pd.concat([all_dat, dat])

    return all_dat


def get_debbink_nonsyn(
        debbink_data_dir = None):
    """
    Get non-synonymous HA substitutions 
    for Debbink et al
    """

    deb_nonsyn = get_debbink(
        debbink_data_dir = debbink_data_dir,
        suffix = "putative.antigenic")
    
    ## find H3 nomenclature and replace
    ## NaN with impossible position value -99999
    deb_nonsyn["H3_pos"] = (
        deb_nonsyn["PDB_4HMG"].str.extract(
            "([0-9]+)").fillna(-99999).astype('int'))
    deb_nonsyn = (
        deb_nonsyn.pipe(compute_antigenic)
    )
    deb_nonsyn["is_aa_substitution"] = True
    deb_nonsyn = deb_nonsyn.rename(
        columns = {"freq.var": "freq_var"})
    excludes = ["season"]
    return_cols = [col for col in column_spec
                   if col in deb_nonsyn.columns
                   and col not in excludes]
    return deb_nonsyn[return_cols]


def get_mccrone(mccrone_isnv_path = None,
                mccrone_nonsyn_path = None):
    """
    Read in non-synonymous
    substitution data from 
    McCrone et al. 2018 
    and do initial cleaning
    """
        
    # make sure keys are strings for good joins
    isnv = pd.read_csv(mccrone_isnv_path,
                       dtype = mc_dtypes)

    nonsyn = pd.read_csv(mccrone_nonsyn_path,
                         dtype = mc_dtypes)

    ## only look at HA nonsynonymous for A/H3N2
    ## and remove mixed infections
    ## this uses a list of mixed specimens per
    ## McCrone et al. published code:
    ## https://github.com/lauringlab/Host_level_IAV_evolution/blob/master/scripts/secondary_analysis/Figures/Figure2.R
    ##
    
    mixed_specimens = ["HS1530", "MH8137", "MH8390"]
    isnv = isnv[(isnv["chr"] == "HA") &
                (isnv["class_factor"] == "Nonsynonymous") &
                (isnv["pcr_result"] == "A/H3N2") &
                (~isnv["SPECID"].isin(mixed_specimens))]
    nonsyn = nonsyn[nonsyn["pcr_result"] == "A/H3N2"]
    nonsyn["is_aa_substitution"] = True

    merge_cols = ["ref", "var", "chr",
                  "class_factor",
                  "freq.var",
                  "AA_pos",
                  "SPECID",
                  "ENROLLID",
                  "onset",
                  "Ref_AA",
                  "Var_AA"]
    
    merge_ids = ["SPECID", "Ref_AA", "Var_AA"]
    additional_cols = ["H3_pos", "is_aa_substitution"]
    result = isnv[merge_cols].merge(
        nonsyn[merge_ids + additional_cols],
        on = merge_ids)
    
    if not result.shape[0] == nonsyn.shape[0]:
        raise RuntimeError("Join failed; check input data.\n"
                           "Path to iSNV data: {}\n"
                           "Path to nonsyn data: {}\n"
                           "nonsyn data shape: {}\n"
                           "joined data shape: {}\n"
                           "".format(mccrone_isnv_path,
                                     mccrone_nonsyn_path,
                                     nonsyn.shape,
                                     result.shape))
    ###########################
    ## clean up amino accids
    ###########################
    ## rename columns to agree with Debbink data spec
    result = result.rename(columns = {"Ref_AA": "aa_consensus",
                                      "Var_AA": "aa_var",
                                      "AA_pos": "aa_pos"})
    clean_AA = lambda x: x[2]
    for col in ["aa_consensus", "aa_var"]:
        result[col] = result[col].apply(clean_AA)
    ## clean up H3 positions, with -99999 as the missing value 
    result["H3_pos"] = result["H3_pos"].fillna(-99999).astype("int")


    result = result.pipe(
        compute_antigenic)
    
    return result
    
    
def get_mccrone_meta(mccrone_meta_path = None):
    """
    Read in metadata for McCrone 
    et al. 2018 and perform initial
    cleaning
    """

    meta = pd.read_csv(mccrone_meta_path,
                       dtype = mc_dtypes)

    ## subset only those who could've shown up in isnv data
    meta_h3 = pd.DataFrame(meta[
        meta.pcr_result.str.contains("H3") &
        meta.snv_qualified &
        meta.sequenced])

    ## create flag for vaccine matching
    meta_h3["vaccination_status"] = (
        meta_h3["vaccination_status"].astype("bool"))
    meta_h3['vaccination_matched'] = (meta_h3['vaccination_status'] &
                                      (meta_h3['season'] != "2014-2015"))

    ## find unique infections (some individual
    ## infections appear as two rows)
    meta_h3 = meta_h3.drop_duplicates(
        ['ENROLLID', 'onset'])
    
    return meta_h3


def make_combined_table(
        debbink_data_dir = None,
        mccrone_isnv_path = None,
        mccrone_nonsyn_path = None,
        mccrone_meta_path = None):
    """
    Read in data and metadata from 
    Debbink et al and McCrone et al,
    and output a clean, complete table 
    for downstream analysis and plotting
    """
    ##############################
    ## handle McCrone et al data
    ##############################

    ## get non-synonymous HA substitutions 
    mc = get_mccrone(
        mccrone_isnv_path = mccrone_isnv_path,
        mccrone_nonsyn_path = mccrone_nonsyn_path)

    ## get metadata 
    mc_meta = get_mccrone_meta(
        mccrone_meta_path = mccrone_meta_path)

    mc_meta["Id"] = mc_meta["ENROLLID"]
    mc["Id"] = mc["ENROLLID"]
    
    ## add on infections without
    ## nonsynonymous HA substitutions 
    mc = mc.sort_values(['is_antigenic_ridge', 
                         'is_antigenic_site',
                         'is_aa_substitution',
                         'Id'],
                        ascending = False
    ).drop_duplicates(
        ["Id", "onset"])

    print("Merging and checking data integrity, McCrone...\n")
    mc = mc_meta.drop_duplicates(
        ["Id"]
    ).merge(
        mc,
        on = "Id",
        how = "left"
    ).rename(
        columns = {"freq.var": "freq_var"}
    ).pipe(
        fill_infection_nas
    ).pipe(
        check_data_integrity)
    
    ## subset just the needed columns
    mc = mc[column_spec]

    ##############################
    ## handle Debbink et al data
    ##############################

    ## get infections with non-synonymous amino acid
    ## substitutions in HA
    deb_nonsyn = get_debbink_nonsyn(
        debbink_data_dir = debbink_data_dir)
    
    ## add on infections with no non-synonymous
    ## substitutions in HA
    deb_all = get_debbink(
        debbink_data_dir = debbink_data_dir,
        suffix = "wg").drop_duplicates(
        "Id")

    ## calculate flag for vaccination matching
    deb_all["vaccination_status"] = deb_all["Vax"]
    deb_all["vaccination_matched"] = (
        (deb_all["vaccination_status"]) &
        (deb_all["season"] != "2004-2005"))

    ## count infections with multiple substitutions
    ## only once, sorting by the most "interesting"
    ## substitution (ridge > site > other)
    all_cols = ["Id",
                "season",
                "vaccination_status",
                "vaccination_matched"]
    print("Merging and checking data integrity, Debbink...\n")
    deb = deb_all[all_cols].merge(
        deb_nonsyn,
        on = "Id",
        how = "left").sort_values(
            ['is_antigenic_ridge', 
             'is_antigenic_site',
             'is_aa_substitution',
             'Id'],
            ascending = False,
    ).drop_duplicates("Id"
    ).pipe(fill_infection_nas
    ).pipe(check_data_integrity)
    
    deb = deb[column_spec]

    result = deb.append(
        mc,
        ignore_index = True,
        verify_integrity = True)

    if not result.shape[0] == deb.shape[0] + mc.shape[0]:
        raise ValueError("combined data frame has different "
                         "number of rows from constituent "
                         "data frames")
    if not ((result.shape[1] == deb.shape[1]) &
            (result.shape[1] == mc.shape[1])):
        raise ValueError("combined data frame has different "
                         "number of columns from constituent "
                         "data frames")

    result["vaccination_class"] = result.apply(
        vaccination_class, axis = 1)
    return result


def main(debbink_data_dir = None,
         mccrone_isnv_path = None,
         mccrone_nonsyn_path = None,
         mccrone_meta_path = None,
         outpath = None):

    table = make_combined_table(
        debbink_data_dir = debbink_data_dir,
        mccrone_isnv_path = mccrone_isnv_path,
        mccrone_nonsyn_path = mccrone_nonsyn_path,
        mccrone_meta_path = mccrone_meta_path)

    table.to_csv(outpath)
    

if __name__ == "__main__":

    if len(sys.argv) < 6:
        print("USAGE ./{} <debbink_data dir> <mccrone isnv data> "
              "<mccrone nonsyn data> <mccrone metadata> "
              "<outpath>".format(sys.argv[0]))
    else:
        main(debbink_data_dir = sys.argv[1],
             mccrone_isnv_path = sys.argv[2],
             mccrone_nonsyn_path = sys.argv[3],
             mccrone_meta_path = sys.argv[4],
             outpath = sys.argv[5])
