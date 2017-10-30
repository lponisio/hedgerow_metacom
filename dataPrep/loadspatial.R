## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/spatial')
source('../../dataPrep/src/misc.R')
library(rgdal)
library(maptools)
library(raster)
save.dir <- '../../data/spatial'

## ************************************************************
## load and prepare data
## ************************************************************
## all hedgerows in landscape
hr.dir <- '../../../hedgerow/spatialData/AllKnownHR'
hr <-
    readOGR(file.path(hr.dir, 'all_hr.shp'),
            'all_hr')

hr.new <- hr[hr@data$Project == "None", "SITE"]
row.names(hr.new) <- hr.new@data$SITE <- paste(hr.new@data$SITE, "new")

colnames(hr.new@data) <- "site"

save(hr, file=file.path(save.dir, "hr.Rdata"))
save(hr.new, file=file.path(save.dir, "hrNew.Rdata"))

## ************************************************************
## leithen's version
## ************************************************************

load('~/Dropbox/hedgerow/data_sets/gis/hedgerows.RData')
line.sites <- dd.h
proj4string(line.sites) <- CRS(proj4string(hr))

save(line.sites,
     file=file.path(save.dir, "sitesLine.Rdata"))

mid.lines <- SpatialLinesMidPoints(line.sites)
sites.pt <- SpatialPointsDataFrame(mid.lines, data=line.sites@data)
proj4string(sites.pt) <- CRS(proj4string(mid.lines))

save(sites.pt,
     file=file.path(save.dir, "sites.Rdata"))

## ************************************************************
## all hedgerows and controls
## ************************************************************

all.sites.lines <- spRbind(hr.new, line.sites)

all.sites.pt <- SpatialLinesMidPoints(all.sites.lines)
all.sites.pt <-  all.sites.pt[!duplicated(all.sites.pt@data$df0),]

## drop the hedgerows that we did not samples that are very far N
## above the others
quartz()
plot(all.sites.pt)
nsites <- nrow(all.sites.pt@data)

## cannot think of anyway but hardcoding this
too.far.N <- all.sites.pt@data$df0[order(coordinates(all.sites.pt)[,2])][nsites:(nsites-14)]

all.sites.pt <-  all.sites.pt[!all.sites.pt@data$df0 %in% too.far.N,]

points(all.sites.pt, col="red")

save(all.sites.lines,
     file=file.path(save.dir, "sitesAllLines.Rdata"))

save(all.sites.pt,
     file=file.path(save.dir, "sitesAll.Rdata"))

## ************************************************************
## yolo landcover
## ************************************************************

landcover <-
    readOGR(path.expand("~/Dropbox/hedgerow/spatialData/yoloLandCover/YoloCounty_RegionalVegetation_July08"),
            "YoloCounty_RegionalVegetation_July08")


landcover <- spTransform(landcover, CRS(proj4string(all.sites.pt)))

non.natural <- c("Water", "Vineyards",
                 "Urban or Built-up",
                 "Truck/Nursery/Berry Crops",
                 "Semiagricultural/Incidental to Agriculture",
                 "Rice",
                 "Pasture",
                 "Levee",
                 "Field Crops",
                 "Deciduous Fruits/Nuts",
                 "Citrus/Subtropical",
                 "Barren - Anthropogenic",
                 "Alkali Sink",
                 "Barren - Gravel and Sand Bars",
                 "Grain/Hay Crops", "Rock Outcrop",
                 "Serpentine Barren",
                 "Giant Reed Series")

kinda.natural <- c("Tamarisk Alliance",
                 "Perennial pepperweed (Lepidium latifolium) Alliance",
                 "Eucalyptus Alliance", "Acacia - Robinia Alliance",
                 "Upland Annual Grasslands & Forbs Formation",
                 "California Annual Grasslands Alliance")


## vegtype <- data.frame(Classification=sort(as.character(unique(landcover@data$VegName))))
## vegtype$Bee.resources <- "yes"
## vegtype$Bee.resources[vegtype$Classification %in% non.natural] <- "no"

## write.table(vegtype,
##             file="~/Dropbox/hedgerow_metacom_saved/occupancy/table/vegtype.txt",
##             sep=" & ", row.names=FALSE)


landcover.nat <- landcover[!landcover@data$VegName %in% non.natural
                          ## & !landcover@data$VegName %in% kinda.natural,]
## include kinda natural components????

landcover.kinda.nat <- landcover[landcover@data$VegName %in% kinda.natural,]

save(landcover.nat, landcover.kinda.nat,
     file=file.path(save.dir, "landcoverNat.Rdata"))



## find out which classes look like they are really ag
f <- function(){
    layout(matrix(1:4, nrow=2))
    for(vegName in unique(landcover.nat$VegName)){
        this.veg <- landcover.nat[landcover.nat$VegName == vegName,]
        plot(this.veg, main=vegName)

    }
}

pdf.f(f, file='../../dataPrep/figures/landcover.pdf')

