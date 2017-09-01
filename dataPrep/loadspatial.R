rm(list=ls())
setwd('analysis/spatial')
## setwd('~/Dropbox/hedgerow_metacom')
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
## all sites?????
## ************************************************************
## site.dir <- "../../../hedgerow/GIS/AllSites"
## sites <-
##   readOGR(file.path(site.dir, 'all_sites.shp'),
##           'all_sites')

## sites.baci <-  sites[sites@data$STUDY == "BACI",]

## save(sites.baci,
##      file=file.path(save.dir, "baciSites.Rdata"))


## save(sites,
##      file=file.path(save.dir, "allSites.Rdata"))


## ************************************************************
## another version...
## ************************************************************
## load('../../data/spatial/sites.Rdata')

## sites <- sites[!duplicated(sites@data$site),]

## save(sites, file='../../data/spatial/sites.Rdata')


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
all.sites.pt <-  all.sites.pt[all.sites.pt@data$df0[!duplicated(all.sites.pt@data$df0)],]

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

non.natural <- c("Water", "Vineyards", "Urban or Built-up",
                 "Truck/Nursery/Berry Crops",
                 "Semiagricultural/Incidental to Agriculture",
                 "Rice",
                 "Pasture", "Levee", "Field Crops",
                 "Deciduous Fruits/Nuts",
                 "Citrus/Subtropical",
                 "Barren - Anthropogenic",
                 "Alkali Sink",
                 "Barren - Gravel and Sand Bars",
                 "Grain/Hay Crops", "Rock Outcrop",
                 "Serpentine Barren",
                 "Tamarisk Alliance", "Eucalyptus Alliance")

landcover.nat <- landcover[!landcover@data$VegName %in% non.natural,]


save(landcover.nat,
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

## make a raster
ill.dim <- ceiling(apply(bbox(landcover.nat), 1, diff)/30)
r <- raster(nrow=ill.dim[1], ncol=ill.dim[2],
            xmn=bbox(landcover.nat)[1,1],
              xmx=bbox(landcover.nat)[1,2],
              ymn=bbox(landcover.nat)[2,1],
              ymx=bbox(landcover.nat)[2,2],
            crs=CRS(proj4string(landcover.nat)))

landcover.r <- rasterize(landcover.nat, r, 'VegCode')

landcover.r <- projectRaster(landcover.r, CRS(proj4string(all.sites.pt)))

save(landcover.r,
     file=file.path(save.dir, "landcoverNat_raster.Rdata"))
