# Asynchrony between virus diversity and antibody selection limits influenza virus evolution
[Dylan H. Morris](https://dylanhmorris.com)(1,\*), [Velislava N. Petrova](2), Fernando W. Rossine(1), [Edyth Parker](https://orcid.org/0000-0001-8312-1446)(3, 4), Bryan T. Grenfell(1, 5), [Richard A. Neher](https://orcid.org/0000-0003-2525-1407)(6), Simon A. Levin(1), and [Colin A. Russell](https://orcid.org/0000-0002-2113-162X)(4,\*)


\* Corresponding authors

1. Department of Ecology \& Evolutionary Biology, Princeton University, Princeton, NJ, USA. 
2. Department of Human Genetics, Wellcome Trust Sanger Institute, Cambridge, UK 
3. Department of Veterinary Medicine, University of Cambridge, Cambridge, UK
4. Department of Medical Microbiology, Academic Medical Center, University of Amsterdam, Amsterdam, NL
5. Fogarty International Center, National Institutes of Health, Bethesda, USA
6. Biozentrum, University of Basel, Basel, CH


## Repository information
This repository accompanies the article "Asynchrony between virus diversity and antibody selection limits influenza virus evolution" (DH Morris et al). It provides code for reproducing numerical simulations and analysis and recreating display figures from the paper.

## License and citation information
If you use the code or data provided here, please make sure to do so in light of the project [license](LICENSE.txt) and please cite our work as below:

- D. H. Morris et al. Asynchrony between virus diversity and antibody selection limits influenza virus evolution. 2020. https://github.com/dylanhmorris/asynchrony-influenza/

Bibtex record:
```
@electronic{morris2020,
    Author = {Dylan H. Morris AND Velislava N. Petrova 
              AND Fernando W. Rossine AND Edyth Parker AND 
              Bryan T. Grenfell AND Richard A. Neher AND
              Simon A. Levin AND Colin A. Russell},
    Title = {Asynchrony between virus diversity and antibody selection limits influenza virus evolution},
    Date = {2020},
    doi = {10.1101/2020.04.27.064915}
}
```

## Article abstract 
Seasonal influenza viruses create a persistent global disease burden by evolving to escape immunity induced by prior infections and vaccinations. New antigenic variants have a substantial selective advantage at the population level, but these variants are rarely selected within-host, even in previously immune individuals. We find that the temporal asynchrony between within-host virus exponential growth and antibody-mediated selection can limit within-host antigenic evolution. If selection for new antigenic variants acts principally at the point of initial virus inoculation, where small virus populations encounter well-matched mucosal antibodies in previously infected individuals, there can exist protection against reinfection that does not regularly produce observable new antigenic variants within individual infected hosts. Our results explain how virus antigenic evolution can be highly selective at the global level but nearly neutral within hosts. They also suggest new avenues for improving influenza control.

## Directories
- ``src``: all code, including data preprocessing, stochastic simulation, numerical analysis, and figure generation:
   - ``src/pysrc``: all Python code, mainly figure generation and data analysis.
   - ``src/gillespie``: Stochastic simulation classes implementing algorithms due to Gillespie and others.
   - ``src/ini``: classes for readining ``.ini`` and ``.mk``-formatted input files. Modifies the open source New BSD licensed [``inih`` library](https://github.com/benhoyt/inih).
   - ``src/random``: pseudorandom number generation classes
   - ``src/testing``: unit tests for C++ code. Uses the Boost-licensed ``catch.hpp`` library.
   - ``src/trees``: code for manipulating and plotting phylogenetic trees.
   - ``src/withinhost``: C++ classes specify within host dynamics models (including transmission chains), as well as C++ classes that hold model parameter sets.
   - ``src/withinhost_runners``: Programs that run stochastic simulations from the models specified in ``src/withinhost``.
- ``dat``: data files in whitespace-separated plain text (``.txt``) and comma-separated values (``.csv``) formats.
- ``out``: output files from simulation, as whitespace-separated plain text (``.txt``) files.
- ``ms``: manuscript files, including main text and extended data figures (as ``.pdf`` files in the ``ms/main/figures/`` directory. The manuscript source ``.tex`` files are not included to respect preprint and embargo policies.
- ``bin``: compiled binaries generated from C++ code, used to run simulations
- ``obj``: object files created during C++ code compilation

## Reproducing analysis

A guide to reproducing the analysis from the paper follows.


### Getting the code
First download this repository. The recommended way is to ``git clone`` it from the command line:

    git clone https://github.com/dylanhmorris/asynchrony-influenza.git

Downloading it manually via Github's download button or from [OSF](https://doi.org/10.17605/OSF.IO/jdsbp/) should also work.

### Dependency installation
The analysis can be auto-run from the project ``Makefile``, but you may need to install some external dependencies first. See the **Dependency installation guide** below for a complete walkthrough. In the first instance, you'll need a working installation of Python 3, a working C++ compiler, and a working installation of Gnu Make or similar. A few external Python packages can then be installed from the command line by typing 

    make depend

from within the project directory.

### Running the analysis

The simplest approach is simply to type ``make`` at the command line, which should produce a full set of figures and simulation output. Be aware, however, that this may take a great deal of time, as some of the simulations, particularly the transmission chains, are computationally costly. Some targets (``pt_chain_[x].txt`` and ``drift_chain_[x].txt``) have been split across multiple files for this reason, so that users with access to multiple cores can parallelize generation by ``make``-ing each in a seperate process. A 4-core parallelized ``make``-ing of the project should take 6--12 hours on a laptop.

If you want to do things piecewise, typing ``make <filename>`` for any of the files listed in the ``out`` or ``ms/main/figures`` directories below should run the steps needed to produce that file.

Some shortcuts are available:

- ``make figs`` produces all figures
- ``make minimal_within_host_results`` produces simulations from the minimal within-host evolutionary model, if they do not already exist
- ``make chain_results`` produces simulations from the transmission chain model, if they do not already exist 
- ``make sensitivity_results`` produces sensitivity analysis simulations from the minimal within-host model, if they do not already exist

- ``make clean`` removes all generated files, leaving only source code (though it does not uninstall packages)
- ``make run_tests`` makes and runs all unit tests for C++ code.

### Examining code

Examining the raw Python and C++ code is the place to start to understand how models have been specified. But note that many parameters are set at runtime rather than hard-coded.

Parameter choices are specified in ``dat/RunParameters.mk``.

## Project structure when complete

Once the full analysis has been run, you should be able to find a full set of figures in ``ms/figures``, simulation output in ``out`` and simulation code binaries in ``bin``.

## Dependency installation guide
You will need a working Python 3 installation with the command line interpreter ``python3`` (macOS and Linux). On mac and Linux, you can check that you have an accessible ``python3`` by typing ``which python3``at the command line and seeing if one is found.

If you do not have an Python 3 installation, you can install it from the [Anaconda project](https://anaconda.org/) or from the command line using a package manager such as [Homebrew](https://brew.sh/) on macOS or ``apt-get`` on Linux. macOS users may also need to install the macOS "command line tools" by typing ``xcode-select --install`` at a command prompt.

To reproduce the phylogenetic tree visualizations in the SI, you will need a working installation of the statistical programmaing language R with the command line interpreter ``Rscript`` (macOS and Linux). On mac and Linux, you can check that you have an accessible ``Rscript`` by typing ``which Rscript``at the command line and seeing if one is found. 

If you do not have an R installation, you can install it from [the R project website](https://www.r-project.org/) or from the command line using a package manager such as [Homebrew](https://brew.sh/) on macOS or ``apt-get`` on Linux. macOS users may also need to install the macOS "command line tools" by typing ``xcode-select --install`` at a command prompt.

Once Python 3 and R are installed, you can automatically install all other dependencies on most systems using ``make``. In the top level project directory, type the following at the command line:

    make depend

Alternatively, you can type the commands ``pip3 install -r python_requirements.txt`` and ``Rscript --vanilla src/trees/install_needed_R_packages.R`` manually. 
