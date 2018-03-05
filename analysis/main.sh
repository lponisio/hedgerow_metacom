#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github

Rscript dataPrep/dataPrep.R

##***************************************************************
## Occupancy analysis
##***************************************************************
## argument: include all interactions or drop those with non-crop
## data (allInt, no_noncrop), non-crop area decay rate (350, 1000,
## 2500), filter over latent states? "filtering, no_filtering"),
## subset the data by site type? (all, hedgerow, control), scale for
## MCMC (# > 0)
Rscript analysis/occupancy/main.R "allInt" "350" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "allInt" "1000" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "allInt" "2500" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "allInt" "350" "filtering" "hedgerow" 1e2
Rscript analysis/occupancy/main.R "allInt" "350" "filtering" "control" 1e2

##***************************************************************
## network roles
##***************************************************************
Rscript analysis/variability/centrality.R
Rscript analysis/variability/plotting/centrality.R
