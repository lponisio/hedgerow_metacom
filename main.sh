#!/usr/bin/env bash

## All paths are relative to the analysis folder within the github
## repo

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github
Rscript dataPrep/loadspatial.R
Rscript dataPrep/dataPrep.R

##***************************************************************
## prep spatial files
##***************************************************************
## number following R scipt represents the number of cores on which to
## run landscape analyses
Rscript analysis/spatial/landscapeHRcalc.R 10

##***************************************************************
## occupancy model
## ***************************************************************
## each model run takes quite a bit of memory and time run each
## combination of renmany habitat decay (first argument, alpha from
## manuscript) and hedgerow decay (second argument).
## note these  models only track top-level parameters by default
Rscript analysis/occupancy/main.R "350" "350" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "1000" "350" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "2500" "350" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "350" "1000" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "1000" "1000" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "1000" "1000" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "350" "2500" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "1000" "2500" "filtering" "all" 1e2
Rscript analysis/occupancy/main.R "2500" "2500" "filtering" "all" 1e2

## drop the super generalists
Rscript analysis/occupancy/main.R "2500" "350" "filtering" "drop.li.ht" 1e2

## visualize results. A warning, sometimes the connection with
## googlemaps doesn't work
Rscript analysis/occupancy/plotting/networks.R

##***************************************************************
## metacommunity network
##***************************************************************
Rscript analysis/networks/spTempMets.R FALSE  "2500" "350"
Rscript analysis/networks/plotting/spTempMets.R FALSE  "2500" "350"

## drop super generalists
Rscript analysis/networks/spTempMets.R TRUE  "2500" "350"
Rscript analysis/networks/plotting/spTempMets.R TRUE  "2500" "350"
