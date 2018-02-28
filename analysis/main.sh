#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github

Rscript ../dataPrep/dataPrep.R

##***************************************************************
## interaction turnover
##***************************************************************
## pollinator partner turnover
Rscript analysis/variability/nulls.R "pols" 99
Rscript analysis/variability/beta-div.R "pols" "abund"
Rscript analysis/variability/beta-div.R "pols" "occ"
Rscript analysis/variability/plotting/beta-div.R "pols" "occ"
Rscript analysis/variability/plotting/beta-div.R "pols" "abund"

## plant partner turnover
Rscript analysis/variability/nulls.R "plants" 99
Rscript analysis/variability/beta-div.R "plants" "abund"
Rscript analysis/variability/beta-div.R "plants" "occ"
Rscript analysis/variability/plotting/beta-div.R "plants" "occ"
Rscript analysis/variability/plotting/beta-div.R "plants" "abund"

##***************************************************************
## network role variability
##***************************************************************
Rscript analysis/variability/centrality.R
Rscript analysis/variability/plotting/centrality.R

##***************************************************************
## phylogenetic signal
##***************************************************************
