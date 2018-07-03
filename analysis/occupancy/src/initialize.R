## ************************************************************
## load needed data sets
## ************************************************************
library(devtools)
library('abind')
library('nimble')
library('R2jags')
source('src/misc.R')
source('src/prep.R')
source('plotting/posteriorPlotting.R')
source('plotting/src/plotting.R')
source('plotting/src/checkChains.R')
source('plotting/src/plotInteractions.R')
source('src/make-matrix.R')
source('src/comparMCMCs_withMonitors.R')
load('../../data/networks/allSpecimens.Rdata')
save.dir <- "../../../hedgerow_metacom_saved/occupancy"

## hedgerow.dir <- "../../../hedgerow/data_sets"
hedgerow.dir <- "../../data"
## spatial data
geo <- read.csv(file.path(hedgerow.dir, 'tables/geography.csv'),
                as.is=TRUE)

## sampling schedule
sr.sched <- read.csv(file.path(hedgerow.dir,
                               'tables/conditions.csv'),
                     as.is=TRUE)
sr.sched$Date <- as.Date(sr.sched$Date)

sr.sched$Site <- geo$Site[match(sr.sched$GeographyFK,
                                geo$GeographyPK)]

## precip data
## precip <- read.csv(file.path(hedgerow.dir,
##                              'misc/precip_between_season_winters.csv'))
## rain <- as.numeric(precip$precip)
## names(rain) <- precip$year

## trait data
all.traits <- read.csv("../../data/traits.csv")
all.traits$BodyLength[!is.na(all.traits$MeanITD)] <-
    all.traits$MeanITD[!is.na(all.traits$MeanITD)]

## HR area of hedgerow in buffers
load('../../data/spatial/HRarea.Rdata')
## HR area weighted by log distance
load('../../data/spatial/HRareaDist.Rdata')
## HR area weighted by gaussian decay
load('../../data/spatial/hrcover_decay.Rdata')
## non-crop area weighted by gaussian decay
load('../../data/spatial/natcover_decay_yolo.Rdata')
## veg data
load('../../data/veg.Rdata')

## raw flower data
load("~/Dropbox/hedgerow/data_sets/traditional/veg-complete.Rdata")
raw.flower.data <- dd


if(length(args) == 0){
    ## allInt, "no_noncrop"
    include.int <- "allInt"
    ## 350, 1000, 2500
    natural.decay  <- "350"
    HR.decay <- "350"
    filtering <- TRUE
    scale <- 1e2
    data.subset <- "all"
}else{
    include.int <- args[1]
    natural.decay <- args[2]
    HR.decay <- args[3]
    filtering <- args[4]
    if(filtering == "filtering"){
        filtering <- TRUE
    } else {
        filtering <- FALSE
    }
    data.subset <- args[5]
    scale <- as.numeric(args[6])
}

if(data.subset == "hedgerow"){
    spec <- spec[spec$SiteStatus == "mature" | spec$SiteStatus ==
                 "maturing",]
}else if(data.subset == "control"){
    spec <- spec[spec$SiteStatus == "control",]
}

sr.sched <- sr.sched[sr.sched$Site %in% unique(spec$Site),]


if(length(args) > 5){
    ncores <- args[7]
}

