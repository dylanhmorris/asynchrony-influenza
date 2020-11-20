#######################################
# filename: Makefile
# author: Dylan Morris <dhmorris@princeton.edu>
#
# description: Makefile for running
# analysis for the Morris et al 
# study of asynchrony between
# diversity and selection in 
# influenza virus evolution
#######################################

#######################################
# commands and expected bash settings 
# 
# (check these match your machine 
# if you are having trouble reproducing
# the analysis)
#######################################

RM = rm -f
CXX = g++-9
CXXFLAGS = -O3 -std=c++14
CXXFLAGSOLD = -O3 -std=c++14
MKDIR = mkdir -p
R_OPTIONS = --vanilla
R_COMMAND := Rscript $(R_OPTIONS)


#################################
# directory structure
##################################

# target locations
OBJDIR := obj
BIN := bin
TEST:= test
MS_SRC := ms
MS_OUT := ms_build
OUTDIR := out

# source code locations
SRCDIR := src
PYTHON := $(SRCDIR)/pysrc
TREE_SRC := $(SRCDIR)/trees

# data locations
DATA := dat
RAW_DATA := $(DATA)/raw
CLEANED_DATA := $(DATA)/cleaned
PHYLO_DATA := $(DATA)/phylo_data

# extensions
SRCEXT := cc
OBJEXT := o
HEADEREXT := h
FIGEXT:= pdf

# standard build locations for binaries
STDLOC = $@

# standard build locations for testing binaries
TESTLOC= $(TEST)/$@

# standard build locations for manuscript
MSLOC = $(MS_OUT)/$@

MS_FIG_DIR = $(MS_SRC)/figures

# subdir names
RUNNERS = runners
TESTING = testing

#################################
# all targets
##################################
.PHONY: all
all: depend run_tests data all_results macros figs


###################################
# external dependency installation
###################################

PY_DEP_FILE = python_requirements.txt
R_DEP_SCRIPT = $(TREE_SRC)/install_needed_R_packages.R

.PHONY: depend
depend:
	pip3 install -r $(PY_DEP_FILE)
	$(R_COMMAND) $(R_DEP_SCRIPT)


#################################
# source variables
#################################


# parameters for various models
RUN_PARAMS = $(DATA)/RunParameters.mk
include $(RUN_PARAMS)

MS_PARAMS = $(DATA)/ManuscriptParameters.mk
include $(MS_PARAMS)


#################################
# data locations
##################################
DEBBINK_DATA_DIR = $(RAW_DATA)/debbink
MCCRONE_ISNV = $(RAW_DATA)/mccrone_isnv.csv
MCCRONE_META = $(RAW_DATA)/mccrone_metadata.csv
MCCRONE_NONSYN = $(RAW_DATA)/mccrone_nonsyn.csv
CLEANED_WH_DATA = $(CLEANED_DATA)/wh_data.csv


#################################
# module objects structure
##################################
GILL_SRC := $(shell find $(SRCDIR)/gillespie -type f -name *.$(SRCEXT))
GILL_OBJ := $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(GILL_SRC:.$(SRCEXT)=.$(OBJEXT)))
GILL_HEADERS := $(shell find $(SRCDIR)/gillespie -type f -name *.$(HEADEREXT))


RNG_SRC := $(shell find $(SRCDIR)/random -type f -name *.$(SRCEXT))
RNG_OBJ:= $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(RNG_SRC:.$(SRCEXT)=.$(OBJEXT)))
RNG_HEADERS := $(shell find $(SRCDIR)/random -type f -name *.$(HEADEREXT))


INI_SRC := $(shell find $(SRCDIR)/ini -type f -name *.$(SRCEXT))
INI_OBJ := $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(INI_SRC:.$(SRCEXT)=.$(OBJEXT)))
INI_HEADERS := $(shell find $(SRCDIR)/ini -type f -name *.$(HEADEREXT))

WH_MODEL_HEADERS := $(shell find $(SRCDIR)/withinhost -type f -name *.$(HEADEREXT))

WITHINHOST_HEADERS = $(GILL_HEADERS) $(RNG_HEADERS) $(INI_HEADERS) $(WH_MODEL_HEADERS)


#################################
# object and binary aliases, for convenience
##################################
WITHINHOST = $(GILL_OBJ) $(RNG_OBJ) $(INI_OBJ)

TEST_NAMES = test_gillespie test_rng test_reader TestMinimalWithinHost TestChainWithinHost TestMinimalParamset
TEST_PATHS = $(addprefix $(TEST)/, $(TEST_NAMES))

PYTEST_NAMES = test_wh_popgen.py
PYTEST_PATHS = $(addprefix $(PYTHON)/, $(PYTEST_NAMES))

FIGURE_BASE_NAMES = plotting_functions.py plotting_style.py


#################################
# data download and cleaning
#################################
DATA_CLEANING_SCRIPT = $(PYTHON)/clean_wh_data.py

$(CLEANED_WH_DATA): $(DATA_CLEANING_SCRIPT) $(DEBBINK_DATA_DIR) $(MCCRONE_ISNV) $(MCCRONE_NONSYN) $(MCCRONE_META)
	$(MKDIR) $(CLEANED_DATA)
	./$^ $@

.PHONY: data

data: $(CLEANED_WH_DATA)

#################################
# transmission chain analyses
#################################

CHAIN_OUTDIR = $(OUTDIR)/transmission_chain_results
CHAIN_RUNNER = TransmissionChain

$(CHAIN_OUTDIR)/pt_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) pt_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/pt_chain_%.txt=%) 0 

$(CHAIN_OUTDIR)/drift_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) drift_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/drift_chain_%.txt=%) 0 

$(CHAIN_OUTDIR)/constant_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) constant_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/constant_chain_%.txt=%) 0 

$(CHAIN_OUTDIR)/better_together_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) better_together_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/better_together_chain_%.txt=%) 0 

$(CHAIN_OUTDIR)/threshold_pt_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) pt_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/threshold_pt_chain_%.txt=%) 1 

$(CHAIN_OUTDIR)/threshold_drift_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) drift_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/threshold_drift_chain_%.txt=%) 1 

$(CHAIN_OUTDIR)/threshold_constant_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) constant_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/threshold_constant_chain_%.txt=%) 1 

$(CHAIN_OUTDIR)/threshold_better_together_chain_%.txt: $(BIN)/$(CHAIN_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(CHAIN_OUTDIR)
	./$< $(RUN_PARAMS) better_together_chain $(SPLIT_N_CHAINS) $(DEFAULT_CHAIN_LENGTH) $@ $(@:$(CHAIN_OUTDIR)/threshold_better_together_chain_%.txt=%) 1


DRIFT_CHAINS = $(addsuffix .txt, $(addprefix drift_chain_, 1 2 3 4))
CONSTANT_CHAINS = $(addsuffix .txt, $(addprefix constant_chain_, 1 2 3 4))
BT_CHAINS = $(addsuffix .txt, $(addprefix better_together_chain_, 1 2 3 4))
PT_CHAINS = $(addsuffix .txt, $(addprefix pt_chain_, 1 2 3 4))

CHAIN_NAMES_STANDARD = $(DRIFT_CHAINS) $(CONSTANT_CHAINS) $(BT_CHAINS) $(PT_CHAINS)
CHAIN_NAMES_THRESHOLD = $(addprefix threshold_, $(CHAIN_NAMES_STANDARD))

CHAIN_NAMES = $(CHAIN_NAMES_STANDARD) $(CHAIN_NAMES_THRESHOLD)

CHAIN_PATHS_STANDARD = $(addprefix $(CHAIN_OUTDIR)/, \
$(CHAIN_NAMES_STANDARD))
CHAIN_PATHS_THRESHOLD = $(addprefix $(CHAIN_OUTDIR)/, \
$(CHAIN_NAMES_THRESHOLD))

CHAIN_PATHS = $(CHAIN_PATHS_STANDARD) $(CHAIN_PATHS_THRESHOLD)

.PHONY: chain_results
chain_results: $(CHAIN_PATHS)

#################################
# effect size distribution analyses
#################################
EFFECT_OUTDIR = $(OUTDIR)/effect_size_results
EFFECT_RUNNER = MutantDistribution


###############################
# within host analyses
###############################
WITHIN_HOST_RESULTS = $(OUTDIR)/within_host_results

MINIMAL_WH_RUNNER = MinimalPSelection
MINIMAL_WH_REPL_RUNNER = MinimalPRepl
MINIMAL_SENSITIVITY_RUNNER = MinimalWhSensitivity

MINIMAL_BN_TARGETS = $(foreach TMP_MIN_MODEL, $(MINIMAL_BN_NAMES), $(addsuffix .txt, $(addprefix $(WITHIN_HOST_RESULTS)/$(TMP_MIN_MODEL)/$(TMP_MIN_MODEL), $(DEFAULT_MWH_BOTTLENECKS))))

MINIMAL_VB_TARGETS = $(foreach TMP_MIN_MODEL, $(MINIMAL_VB_NAMES), $(addsuffix .txt, $(addprefix $(WITHIN_HOST_RESULTS)/$(TMP_MIN_MODEL)/$(TMP_MIN_MODEL), $(DEFAULT_MWH_VB_RATIOS))))

MINIMAL_TARGETS = $(MINIMAL_BN_TARGETS) $(MINIMAL_VB_TARGETS)

MINIMAL_TEST_TARGETS = $(addsuffix .txt, $(addprefix $(WITHIN_HOST_RESULTS)/$(MINIMAL_TEST_NAME)/$(MINIMAL_TEST_NAME), $(DEFAULT_MWH_VB_RATIOS))))

.PHONY: minimal_within_host_results minimal_test_results
minimal_within_host_results: $(MINIMAL_TARGETS)
minimal_test_results: $(MINIMAL_TEST_TARGETS)

$(WITHIN_HOST_RESULTS)/$(MINIMAL_STERILIZING_NAME)/$(MINIMAL_STERILIZING_NAME)%.txt: $(BIN)/$(MINIMAL_WH_RUNNER)
	$(MKDIR) $(dir $@)
	$(eval TMP_BOTTLENECK:= $(patsubst $(MINIMAL_STERILIZING_NAME)%, %, $(basename $(notdir $@))))
	./$(BIN)/$(MINIMAL_WH_RUNNER) $(RUN_PARAMS) $(MINIMAL_STERILIZING_NAME) $(MINIMAL_STERILIZING_NITER) $(TMP_BOTTLENECK) $(basename $@) 0


# define a make rule for all other minimal model results
$(WITHIN_HOST_RESULTS)/%.txt: $(BIN)/$(MINIMAL_WH_RUNNER)
	$(MKDIR) $(dir $@)
	$(eval TMP_MODEL_NAME := $(patsubst %/,%,$(dir $(@:$(WITHIN_HOST_RESULTS)/%.txt=%))))
	$(eval TMP_VB_RATIO := $(patsubst $(TMP_MODEL_NAME)%, %, $(basename $(notdir $@))))
	./$(BIN)/$(MINIMAL_WH_RUNNER) $(RUN_PARAMS) $(TMP_MODEL_NAME) $(MINIMAL_NITER) $(TMP_VB_RATIO) $(basename $@) 0


###############################
# sensitivity analyses
###############################
SENSITIVITY_RESULTS = $(OUTDIR)/sensitivity_analysis_results

SENSITIVITY_WH_RUNNER = MinimalWhSensitivity

SENSITIVITY_TARGETS = $(foreach TMP_MODEL, $(SENSITIVITY_NAMES), $(addsuffix .txt, $(addprefix $(SENSITIVITY_RESULTS)/$(TMP_MODEL)/$(TMP_MODEL), $(SENSITIVITY_BOTTLENECKS))))


$(SENSITIVITY_RESULTS)/%.txt: $(BIN)/$(SENSITIVITY_WH_RUNNER) #$(RUN_PARAMS)
	$(MKDIR) $(dir $@)
	$(eval TMP_MODEL_NAME := $(patsubst %/,%,$(dir $(@:$(SENSITIVITY_RESULTS)/%.txt=%))))
	$(eval TMP_BOTTLENECK := $(patsubst $(TMP_MODEL_NAME)%, %, $(basename $(notdir $@))))
	./$(BIN)/$(SENSITIVITY_WH_RUNNER) $(RUN_PARAMS) $(TMP_MODEL_NAME) $(SENSITIVITY_NITER) $(SENSITIVITY_N_PARAMSETS) $(TMP_BOTTLENECK) $(basename $@) 0

.PHONY: sensitivity_results
sensitivity_results: $(SENSITIVITY_TARGETS)


###############################
# replication selection analysis
###############################
REPL_RESULTS = $(OUTDIR)/replication_selection_results
REPL_PREFIX = p_repl

## paths of replication selection
## simulations (cartesian product
## of k values and t_M values)
REPL_PATHS := $(foreach K, $(REPL_KS), $(foreach T_M, $(REPL_T_MS), $(REPL_RESULTS)/$(REPL_PREFIX)-$(K)-$(T_M).txt)) 

$(REPL_RESULTS)/%.txt: $(BIN)/$(MINIMAL_WH_REPL_RUNNER)
	$(MKDIR) $(dir $@)
	$(eval TMP_K_T_M := $(subst -, , $(patsubst $(REPL_RESULTS)/$(REPL_PREFIX)-%.txt, %, $@)))
	$(eval TMP_K := $(word 1, $(TMP_K_T_M)))
	$(eval TMP_T_M := $(word 2, $(TMP_K_T_M)))
	./$< $(RUN_PARAMS) $(MINIMAL_VISVAR_NAME) $(REPL_NITER) $(REPL_T_FINAL) $(DEFAULT_BOTTLENECK) $(TMP_K) $(TMP_T_M) $@ 0

.PHONY: repl_results
repl_results: $(REPL_PATHS)


###############################
# all analyses
###############################
.PHONY: all_results

all_results: $(REPL_PATHS) $(SENSITIVITY_TARGETS) $(MINIMAL_TARGETS) $(CHAIN_PATHS)


############
# Figures
############
TREE_FIGURE_NAMES = figure-polyphyly-H1.pdf figure-polyphyly-H3-2008-2011.pdf figure-polyphyly-H3-2012-2014.pdf

MAIN_FIG_NAMES := figure-wh-dynamics-summary.pdf figure-inoculation-summary.pdf figure-mutant-dists.pdf figure-population-level-summary.pdf figure-minimal-immune-compromised.pdf figure-escape-tradeoff.pdf figure-wh-sterilizing-timecourse.pdf figure-sensitivity.pdf figure-emergence-time-cdf.pdf figure-prepl-analytic-vs-sims.pdf figure-sensitivity-scatter-fixed.pdf figure-sensitivity-scatter-visvar.pdf figure-replicator-equation.pdf figure-mutant-dists-threshold.pdf figure-relative-contribution.pdf figure-smith-antigenic-map.pdf $(TREE_FIGURE_NAMES)

MS_FIGS = $(addprefix $(MS_FIG_DIR)/, $(MAIN_FIG_NAMES))

FIGURE_SCRIPTS_DIR = $(SRCDIR)/pysrc
FIG_SCRIPTS = $(SRCDIR)/pysrc

FIG_POP_SCRIPTS := figure_analytic_escape_prob.py figure_epimodel_quantiles.py
FIG_POP_DEPS := $(patsubst %.py, $(FIG_SCRIPTS)/%.py, $(FIG_POP_SCRIPTS))

FIG_BASE_SCRIPTS = plotting_style.py plotting_functions.py
FIGURE_BASE := $(patsubst %.py, $(FIG_SCRIPTS)/%.py, $(FIG_BASE_SCRIPTS))

.PHONY: figs
figs: $(MS_FIGS)

WH_DYNAMICS_DEPS = $(FIGURE_SCRIPTS_DIR)/figure_wh_dynamics_summary.py $(FIGURE_BASE) $(MINIMAL_VB_TARGETS) $(CLEANED_WH_DATA) $(FIG_SCRIPTS)/figure_wh_ngs.py
$(MS_FIG_DIR)/figure-wh-dynamics-summary.pdf: $(WH_DYNAMICS_DEPS)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(WITHIN_HOST_RESULTS)/$(MINIMAL_FIXED_NAME) $(WITHIN_HOST_RESULTS)/$(MINIMAL_VISVAR_NAME) $(CLEANED_WH_DATA) $(RUN_PARAMS) $@ $(EXAMPLE_TIMECOURSE_BOTTLENECK) $(DEFAULT_BOTTLENECK)

INOC_SUMMARY_ARGS = $(DATA)/bottleneck-schematic.png $(WITHIN_HOST_RESULTS)/$(MINIMAL_VISVAR_NAME) $(RUN_PARAMS)

INOC_SUMMARY_EXTRA_DEPS = $(FIGURE_BASE) $(FIGURE_SCRIPTS_DIR)/figure_distribution_shift.py $(FIGURE_SCRIPTS_DIR)/figure_optimal_selector.py $(MINIMAL_VB_TARGETS)

$(MS_FIG_DIR)/figure-inoculation-summary.pdf: $(FIGURE_SCRIPTS_DIR)/figure_inoculation_summary.py $(INOC_SUMMARY_ARGS) $(INOC_SUMMARY_EXTRA_DEPS)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(INOC_SUMMARY_ARGS) $@ 

$(MS_FIG_DIR)/figure-prepl-analytic-vs-sims.pdf: $(FIGURE_SCRIPTS_DIR)/figure_prepl_analytic_vs_sims.py $(REPL_RESULTS) $(REPL_PATHS) $(FIGURE_BASE)
	$(MKDIR) $(MS_FIG_DIR)
	./$< $(REPL_RESULTS) $@

$(MS_FIG_DIR)/figure-mutant-dists.pdf: $(FIGURE_SCRIPTS_DIR)/figure_mutant_dists.py $(FIGURE_BASE) $(CHAIN_PATHS_STANDARD) 
	@mkdir -p $(MS_FIG_DIR)
	./$< $(CHAIN_PATHS_STANDARD) $(RUN_PARAMS) $@

$(MS_FIG_DIR)/figure-mutant-dists-threshold.pdf: $(FIGURE_SCRIPTS_DIR)/figure_mutant_dists.py $(FIGURE_BASE) $(RUN_PARAMS) $(CHAIN_PATHS_THRESHOLD) 
	@mkdir -p $(MS_FIG_DIR)
	./$< $(CHAIN_PATHS_THRESHOLD) $(RUN_PARAMS) $@


$(MS_FIG_DIR)/figure-population-level-summary.pdf: $(FIGURE_SCRIPTS_DIR)/figure_population_level_summary.py $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(RUN_PARAMS) $@

$(MS_FIG_DIR)/figure-minimal-immune-compromised.pdf: $(FIGURE_SCRIPTS_DIR)/figure_minimal_immune_compromised.py $(RUN_PARAMS) $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(DEFAULT_BOTTLENECK) $(RUN_PARAMS) $@

$(MS_FIG_DIR)/figure-escape-tradeoff.pdf: $(FIGURE_SCRIPTS_DIR)/figure_escape_tradeoff.py $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(RUN_PARAMS) $@

$(MS_FIG_DIR)/figure-wh-sterilizing-timecourse.pdf: $(FIGURE_SCRIPTS_DIR)/figure_wh_sterilizing_timecourse.py $(MINIMAL_BN_TARGETS) $(RUN_PARAMS) $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(WITHIN_HOST_RESULTS)/$(MINIMAL_STERILIZING_NAME) $(RUN_PARAMS) $@ 

$(MS_FIG_DIR)/figure-sensitivity.pdf: $(FIGURE_SCRIPTS_DIR)/figure_sensitivity.py $(SENSITIVITY_TARGETS) $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(SENSITIVITY_RESULTS) $@ 

$(MS_FIG_DIR)/figure-sensitivity-scatter-visvar.pdf: $(FIGURE_SCRIPTS_DIR)/figure_sensitivity_scatter.py $(SENSITIVITY_TARGETS) $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(SENSITIVITY_RESULTS)/minimal_visvar $(RUN_PARAMS) $@

$(MS_FIG_DIR)/figure-sensitivity-scatter-fixed.pdf: $(FIGURE_SCRIPTS_DIR)/figure_sensitivity_scatter.py $(SENSITIVITY_TARGETS) $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(SENSITIVITY_RESULTS)/minimal_fixed $(RUN_PARAMS) $@

$(MS_FIG_DIR)/figure-smith-antigenic-map.pdf: $(DATA)/smith-antigenic-map.pdf
	@mkdir -p $(MS_FIG_DIR)
	cp $(DATA)/smith-antigenic-map.pdf $(MS_FIG_DIR)/figure-smith-antigenic-map.pdf


$(MS_FIG_DIR)/figure-emergence-time-cdf.pdf: $(FIGURE_SCRIPTS_DIR)/figure_emergence_time_cdf.py $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $(WITHIN_HOST_RESULTS)/$(MINIMAL_FIXED_NAME) $(WITHIN_HOST_RESULTS)/$(MINIMAL_VISVAR_NAME) $(RUN_PARAMS) $@

$(MS_FIG_DIR)/figure-replicator-equation.pdf: $(FIGURE_SCRIPTS_DIR)/figure_replicator_equation.py $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $@

$(MS_FIG_DIR)/figure-relative-contribution.pdf: $(FIGURE_SCRIPTS_DIR)/figure_relative_contribution.py $(FIGURE_BASE)
	@mkdir -p $(MS_FIG_DIR)
	./$< $@ 


# polyphyletic emergence figures
############

TREE_DATA_NAMES = H1N1.nwk H3N2_2008_2011.nwk H3N2_2012_2014.nwk
TREE_POSITION_DATA_NAMES = H1_position_140.txt H3_08_position_158_189.txt H3_12_position_159.txt

TREE_PLOTTER = $(TREE_SRC)/plot_trees.R
TREE_DATA_PATHS =  $(addprefix $(PHYLO_DATA)/, $(TREE_DATA_NAMES))
TREE_POSITION_DATA =  $(addprefix $(PHYLO_DATA)/, $(TREE_POSITION_DATA_NAMES))

TREE_FIGURE_PATHS =  $(addprefix $(MS_FIG_DIR)/, $(TREE_FIGURE_NAMES))

$(MS_FIG_DIR)/figure-polyphyly-%.pdf: $(TREE_PLOTTER) $(TREE_DATA_PATHS) $(TREE_POSITION_DATA)
	@mkdir -p $(MS_FIG_DIR)
	$(R_COMMAND) $^ $(TREE_FIGURE_PATHS)
	@$(RM) Rplots.pdf # clean up Rscript extras

###########################
# dynamic parameter macros
# (used to autopopulate
# quantities cited in the
# manuscript via LaTeX
# macros)
###########################
MACRO_GENERATOR_NAME = make_parameter_macros.py
MACRO_SCRIPTS_DIR = $(PYTHON)/dynamic_results
MACRO_GENERATOR_PATH = $(MACRO_SCRIPTS_DIR)/$(MACRO_GENERATOR_NAME)
MAIN_MACRO_DIR = $(MS_SRC)

MAIN_PARAMETER_MACROS := $(MAIN_MACRO_DIR)/parameters.sty
PARAMETER_MACRO_TEMPLATE = $(DATA)/parameters.sty.tmpl

$(MAIN_PARAMETER_MACROS): $(MACRO_GENERATOR_PATH) $(RUN_PARAMS) $(PARAMETER_MACRO_TEMPLATE)
	$(MKDIR) $(MAIN_MACRO_DIR)
	./$(MACRO_GENERATOR_PATH) $(RUN_PARAMS) $(PARAMETER_MACRO_TEMPLATE) $@

.PHONY: macros
macros: $(MAIN_PARAMETER_MACROS)


#############################
# rules for binaries
#############################

$(BIN)/$(MINIMAL_WH_RUNNER): $(OBJDIR)/$(RUNNERS)/$(MINIMAL_WH_RUNNER).o $(WITHINHOST)
	@$(MKDIR) $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $(STDLOC)

$(BIN)/$(MINIMAL_SENSITIVITY_RUNNER): $(OBJDIR)/$(RUNNERS)/$(MINIMAL_SENSITIVITY_RUNNER).o $(WITHINHOST)
	@$(MKDIR) $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $(STDLOC)

$(BIN)/$(MINIMAL_WH_REPL_RUNNER): $(OBJDIR)/$(RUNNERS)/$(MINIMAL_WH_REPL_RUNNER).o $(WITHINHOST)
	@$(MKDIR) $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $(STDLOC)


$(BIN)/MinimalSterilizing: $(OBJDIR)/$(RUNNERS)/MinimalSterilizing.o $(WITHINHOST)
	@$(MKDIR) $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $(STDLOC)


$(BIN)/TransmissionChain: $(OBJDIR)/$(RUNNERS)/TransmissionChain.o $(WITHINHOST) $(SRCDIR)/$(RUNNERS)/ChainModel.h
	@$(MKDIR) $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $(STDLOC)


$(BIN)/MutantDistribution: $(OBJDIR)/$(RUNNERS)/MutantDistribution.o $(WITHINHOST)
	@$(MKDIR) $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $(STDLOC)

#############################
# rules for test binaries
#############################

$(TEST)/test_gillespie: $(OBJDIR)/$(TESTING)/test_gillespie.o $(RNG_OBJ) $(GILL_OBJ)
	@$(MKDIR) $(TEST)
	$(CXX) $(CXXFLAGSOLD) $^ -o $@

$(TEST)/test_rng: $(OBJDIR)/$(TESTING)/test_rng.o $(RNG_OBJ)
	@$(MKDIR) $(TEST)
	$(CXX) $(CXXFLAGSOLD) $^ -o $@

$(TEST)/test_reader: $(OBJDIR)/$(TESTING)/test_reader.o $(INI_OBJ)
	@$(MKDIR) $(TEST)
	$(CXX) $(CXXFLAGSOLD) $^ -o $@

$(TEST)/TestMinimalWithinHost: $(OBJDIR)/$(TESTING)/TestMinimalWithinHost.o $(WITHINHOST)
	@$(MKDIR) $(TEST)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST)/TestChainWithinHost: $(OBJDIR)/$(TESTING)/TestChainWithinHost.o $(WITHINHOST)
	@$(MKDIR) $(TEST)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(TEST)/TestMinimalParamset: $(OBJDIR)/$(TESTING)/TestMinimalParamset.o $(WITHINHOST)
	@$(MKDIR) $(TEST)
	$(CXX) $(CXXFLAGS) $^ -o $@


#############################
# dependency rules
#############################

$(OBJDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@


###################################
# convenience group rules / .phony
###################################
DIRS := . $(shell find $(SRCDIR) -type d)
GARBAGE_PATTERNS := *~
GARBAGE := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(GARBAGE_PATTERNS)))

.PHONY: test run_tests clean garbage_collect clobber

# clean up object code and executables
clean:
	$(RM) -r $(BIN)
	$(RM) -r $(TEST)
	$(RM) -r $(OBJDIR)
	$(RM) -r $(MSOUT)

# clean up emacs temp files, including in subdirs
garbage_collect:
	@$(RM) $(GARBAGE)

clobber: clean garbage_collect

# conduct unit testing
test: $(TEST_PATHS) $(PYTEST_PATHS)

run_tests: test
	for prog in $(TEST_NAMES); do ./$(TEST)/$$prog; done
	for prog in $(PYTEST_NAMES); do ./$(PYTHON)/$$prog; done


