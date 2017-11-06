#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github

Rscript dataPrep/dataPrep.R
Rscript dataPrep/loadspatial.R

##***************************************************************
## prep spatial files
##***************************************************************
Rscript analysis/spatial/landscapeHRcalc.R

##***************************************************************
## occupancy model
##***************************************************************
Rscript analysis/occupancy/main.R

## visualize results
Rscript analysis/occupancy/plotting/posteriorPlotting.R
Rscript analysis/occupancy/plotting/networks.R
