## ************************************************************
## load needed data sets
## ************************************************************
library('abind')
library('nimble')
source('src/misc.R')
source('src/prep.R')
source('src/setup.R')
source('plotting/posteriorPlotting.R')
source('plotting/src/plotting.R')
source('plotting/src/checkChains.R')
source('plotting/src/plotInteractions.R')
source('plotting/src/plotResults.R')
source('src/make-matrix.R')
source('src/tests.R')
load('../../data/networks/allSpecimens.Rdata')

save.dir <- "../../../hedgerow_metacom_saved/occupancy"
checkDirExists(save.dir)
checkDirExists(file.path(save.dir, "runs"))

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

## trait data
all.traits <- read.csv("../../data/traits.csv")
all.traits$MeanITD[!is.na(all.traits$MeanITD)] <-
    all.traits$MeanITD[!is.na(all.traits$MeanITD)]

## HR area weighted by gaussian decay
load('../../data/spatial/hrcover_decay.Rdata')
## ## non-crop area weighted by gaussian decay
load('../../data/spatial/natcover_decay_yolo.Rdata')
## veg data
load('../../data/veg.Rdata')

if(length(args) == 0){
    natural.decay  <- "2500"
    HR.decay <- "350"
    filtering <- TRUE
    scale <- 1e2
    data.subset <- "all"
}else{
    natural.decay <- args[1]
    HR.decay <- args[2]
    filtering <- args[3]
    if(filtering == "filtering"){
        filtering <- TRUE
    } else {
        filtering <- FALSE
    }
    data.subset <- args[4]
    scale <- as.numeric(args[5])
}

if(data.subset == "hedgerow"){
    spec <- spec[spec$SiteStatus == "mature" | spec$SiteStatus ==
                 "maturing",]
}
if(data.subset == "control"){
    spec <- spec[spec$SiteStatus == "control",]
}
if(data.subset == "drop.li.ht"){
    spec <- spec[!spec$GenusSpecies %in%
                 c("Lasioglossum (Dialictus) incompletum",
                   "Halictus tripartitus"), ]
}
sr.sched <- sr.sched[sr.sched$Site %in% unique(spec$Site),]

if(filtering){
    source(sprintf('src/models/filter.R'))
} else {
    source(sprintf('src/models/latent.R'))
}
