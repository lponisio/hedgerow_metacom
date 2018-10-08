## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/spatial')
source('../../dataPrep/src/misc.R')
source('../../dataPrep/src/plotLandCovers.R')
library(rgdal)
library(maptools)
library(raster)
library(dismo)
library(viridis)
library(rgeos)
library(spectralGP)
library(sp)
save.dir <- '../../data/spatial'
load('../../data/networks/allSpecimens.Rdata')

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
dists <- rdist.earth(coordinates(all.sites.pt), miles=FALSE)

## any sites within 10 m of each other?
apply(dists < 0.01, 1, sum)
## only the diagonal...


## drop the hedgerows that we did not samples that are very far N
## above the others
quartz()
plot(all.sites.pt)
nsites <- nrow(all.sites.pt@data)

## cannot think of anyway but hardcoding this
too.far.N <- all.sites.pt@data$df0[order(coordinates(all.sites.pt)[,2])][nsites:(nsites-14)]

all.sites.pt <-  all.sites.pt[!all.sites.pt@data$df0 %in% too.far.N,]
all.sites.lines <-  all.sites.lines[!all.sites.lines@data$site %in% too.far.N,]

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
                   "California Annual Grasslands Alliance",
                   "Upland Annual Grasslands & Forbs Formation")

## make table for ms
vegtype <- data.frame(Classification=sort(as.character(unique(landcover@data$VegName))))
vegtype$Bee.resources <- "yes"
vegtype$Bee.resources[vegtype$Classification %in% non.natural] <- "no"

write.table(vegtype,
            file="~/Dropbox/hedgerow_metacom_saved/occupancy/table/vegtype.txt",
            sep=" & ", row.names=FALSE)

landcover.nat <- landcover[!landcover@data$VegName %in% non.natural,]

new.hr <- all.sites.pt[grep("new", all.sites.pt@data$df0),]
surveyed.hr <- all.sites.pt[all.sites.pt@data$df0 %in%
                            unique(spec$Site[spec$SiteStatus %in% c("maturing", "mature")]),]

## find out which classes look like they are really ag
f <- function(){
    ## layout(matrix(1:4, nrow=2))
    for(vegName in unique(landcover.nat$VegName)){
        this.veg <- landcover.nat[landcover.nat$VegName == vegName,]
        plot(this.veg, main=vegName)
        points(new.hr, col="red")
        points(surveyed.hr, col="dodgerblue")


    }
}

pdf.f(f, file="../../../hedgerow_metacom_saved/map/landcovers.pdf")

## merge polygones into single "natural" layer
landcover.nat@data$type <- "natural"
landcover.nat <- unionSpatialPolygons(landcover.nat,
                                      landcover.nat@data$type)

dats <- data.frame(type="natural")
rownames(dats) <- "natural"
landcover.nat <- SpatialPolygonsDataFrame(landcover.nat,
                                          data=dats)
save(landcover.nat,
     file=file.path(save.dir, "landcoverNat.Rdata"))

## ************************************************************
## plotting
## ************************************************************
### mered polygons by adj and type
landcover.nat <-
    readOGR(path.expand("~/Dropbox/hedgerow/spatialData/yoloLandCover/landcover_Dissolve_byVegName"),
            "landcover_Dissolve_byVegName_notAdjacent")

pdf.f(plotLandscapeCovers,
      file="../../../hedgerow_metacom_saved/map/covermap.pdf",
      height=10, width=8)

save(landcover.nat,
     file=file.path(save.dir, "landcoverNatPatches.Rdata"))

## ************************************************************
## calculating remant "patch" sizes
## ************************************************************
landcover.patch <-
    readOGR(path.expand("~/Dropbox/hedgerow/spatialData/yoloLandCover/landcover_merged_patches"),
            "landcover_dissolve_nearbyPolys")

landcover.patch@data$area <- rgeos::gArea(landcover.patch, byid=TRUE)

plotPatchHist <- function(){
    par(oma=c(3,3,1,1))
    barplot(table(round(landcover.patch@data$area, -4)),
            las=2)
    mtext("Frequency", 2, line=4, cex=1.2)
    mtext("Area (m^2)", 1, line=5,cex=1.2)
}

pdf.f(plotPatchHist,
      file="../../../hedgerow_metacom_saved/map/patchHist.pdf",
      height=4, width=8)

## ************************************************************
## looking at each type of vegetation on the landscape
## ************************************************************

## landcover <- spTransform(landcover,
##                              CRS("+init=epsg:4326"))

## landcover <- spTransform(landcover,
##                              sys@crs)


## cols <- add.alpha(viridis(length(unique(landcover@data$VegName))), 0.5)
## names(cols) <- unique(landcover@data$VegName)

## for(veg.type in unique(landcover@data$VegName)){
##     print(veg.type)
##     plot(landcover[landcover@data$VegName == veg.type,],
##          col=cols[veg.type], add=TRUE)
## }
