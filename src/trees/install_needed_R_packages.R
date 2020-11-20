#!/usr/bin/env Rscript

#####################################
## name: install_needed_R_packages.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## installs needed R packages for
## reproducing phylogenetic tree
## figures for Morris et al
## paper on asynchrony between
## diversity and selection in ifluenza
####################################

install_if_absent <- function(package_name){
    if (!suppressPackageStartupMessages(
             require(package_name, character.only=TRUE))){
      install.packages(pkgs=package_name,
                       repos="http://cloud.r-project.org")
  }
  else
      cat(sprintf("Package %s already installed\n", package_name))
}

needed_packages <- c(
    "ape",
    "dplyr",
    "ggplot2",
    "cowplot",
    "colorspace",
    "BiocManager")

for (package in needed_packages)
    install_if_absent(package)

## ggtree requires special installation via Bioconductor
if (!suppressPackageStartupMessages(
         require("ggtree", character.only=TRUE))){
    BiocManager::install("ggtree")
} else {cat("Package ggtree already installed\n")}
