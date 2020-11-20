# Data 

## Meta-analysis
The ``raw`` folder contains unprocessed data for the meta-analysis, based on two papers:

- Debbink et al. 2017. Vaccination has minimal impact on the intrahost diversity of H3N2 influenza viruses. PLoS Pathogens. https://doi.org/10.1371/journal.ppat.1006194

- McCrone et al. 2018. Stochastic processes constrain the within and between host evolution of influenza virus. eLife. https://doi.org/10.7554/eLife.35962.001

- ``debbink``: Data files by year from Debbink et al. downloaded from the authors' [published data on Github](https://github.com/lauringlab/Fluvacs_paper/tree/master/results)
- ``mccrone_metadata.csv``: downloaded from the original authors' [published data on Github](https://raw.githubusercontent.com/lauringlab/Host_level_IAV_evolution/master/data/meta_snv_qual.csv); original filename: ``meta_snv_qual.csv``
- ``mccrone_isnv.csv``: downloaded from authors' [published data on Github](https://raw.githubusercontent.com/lauringlab/Host_level_IAV_evolution/master/data/processed/secondary/qual.snv.csv); original filename: ``qual.snv.csv``
- ``mccrone_nonsyn.csv``: downloaded from [data published alongside article in eLife](https://doi.org/10.7554/eLife.35962.015); original filename: ``elife-35962-fig2-data2-v3.csv``.

This data is then processed further to form the file ``cleaned/clean_wh_data.csv``

## Phylogenetic analysis
The ``pylo_data`` folder contains data for the phylogenetic analysis comes from original phylogenies built from [GISAID](https://www.gisaid.org/) EpiFlu(tm) surveillance data, per phylogenetic methods described in the SI. 
The ``gisaid`` folder contains accession numbers for sequences used, with credit to contributing labs.

## Run Parameters
``RunParameters.mk`` contains all parameter settings for stochastic simulations and numerically computed models. The parameters for a particular model are preceded by its model name, which will typically correspond to a file and/or folder in the ``out`` directory. Model names overload defaults. If we are simulating from the model ``minimal_visvar``, for example, the parameter $t_M$ will default to ``DEFAULT_T_M`` unless we set ``MINIMAL_VISVAR_T_M``.

## Other files
Two other files are unlikely to be of use in reproducing the analysis, but are provided for your reference:
- ``ManuscriptParameters.mk`` contains parameters that are not for modeling but simply for display (e.g. the fineness of the heatmap in Fig. 1 of the main text.
- ``parameters.sty.tmpl`` is a Jinja2 template that is filled by the ``src/pysrc/dynamic_results/make_parameter_macros.py`` script (type ``make macros`` to run it). This enables model parameters and output values to be referred to directly in the $\LaTeX$ source for the manuscript. 
