## ************************************************************
## load needed data sets
## ************************************************************
source('src/make-matrix.R')
source('src/comparMCMCs_withMonitors.R')
load('../../data/networks/allSpecimens.Rdata')
save.dir <- "../../../hedgerow_metacom_saved/occupancy"

hedgerow.dir <- "../../../hedgerow/data_sets"
## spatial data
geo <- read.csv(file.path(hedgerow.dir, 'tables/geography.csv'),
                as.is=TRUE)

## sampling schedule
sr.sched <- read.csv(file.path(hedgerow.dir, 'tables/conditions.csv'),
                as.is=TRUE)
sr.sched$Date <- as.Date(sr.sched$Date)

sr.sched$Site <- geo$Site[match(sr.sched$GeographyFK,
                                geo$GeographyPK)]

## precip data
precip <- read.csv(file.path(hedgerow.dir,
                             'misc/precip_between_season_winters.csv'))
rain <- as.numeric(precip$precip)
names(rain) <- precip$year

## trait data
all.traits <- read.csv("../../data/traits.csv")
all.traits$Nestuse <- paste(all.traits$NestLoc, all.traits$Excavate)
all.traits$Nestuse[all.traits$Nestuse == "NA NA"] <- NA
all.traits$Nestuse[all.traits$Nestuse == "NA rent"] <- NA
all.traits$BodyLength[!is.na(all.traits$MeanITD)] <-
    all.traits$MeanITD[!is.na(all.traits$MeanITD)]

## area of hedgerow in buffers
load('../../data/spatial/HRarea.Rdata')
## area weighted by log distance
load('../../data/spatial/HRareaDist.Rdata')
## area of different landcovers
load('../../data/spatial/natcover_decay_yolo.Rdata')
load('../../data/spatial/kinda_natcover_decay_yolo.Rdata')
## veg data
load('../../data/veg.Rdata')


