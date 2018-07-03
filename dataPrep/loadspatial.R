## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/spatial')
source('../../dataPrep/src/misc.R')
library(rgdal)
library(maptools)
library(raster)
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
## landcover.nat@data$type <- "natural"
## landcover.nat <- unionSpatialPolygons(landcover.nat,
##                                       landcover.nat@data$type)

## dats <- data.frame(type="natural")
## rownames(dats) <- "natural"
## landcover.nat <- SpatialPolygonsDataFrame(landcover.nat,
##                                           data=dats)

save(landcover.nat,
     file=file.path(save.dir, "landcoverNat.Rdata"))




## plotting
library(dismo)
library(viridis)
plotLandscapeCovers <- function(){
    ## prep map and site data
    all.sites.pt <- spTransform(all.sites.pt,
                                CRS("+init=epsg:4326"))

    landcover.nat <- spTransform(landcover.nat,
                                 CRS("+init=epsg:4326"))

    bbox.sites <- bbox(all.sites.pt)
    bbox.nat <- bbox(landcover.nat)


    bbox.all <- matrix(c(bbox.sites[1,1],
                         bbox.sites[2,1],
                         bbox.sites[1,2],
                         bbox.sites[2,2]),
                       ncol=2)


    sys <- gmap(x=bbox.all,
                scale=2,
                type="satellite", zoom=10)

    all.sites.pt <- spTransform(all.sites.pt,
                                sys@crs)

    landcover.nat <- spTransform(landcover.nat,
                                 sys@crs)


    new.hr <- all.sites.pt[grepl("new", all.sites.pt@data$df0),]
    surveyed.hr <- all.sites.pt[all.sites.pt@data$df0 %in%
                                unique(spec$Site[spec$SiteStatus %in%
                                                 c("maturing", "mature")]),]
    surveyed.control <- all.sites.pt[all.sites.pt@data$df0 %in%
                          unique(spec$Site[spec$SiteStatus =="control"]),]

    dims <- bbox(sys)
    plot(sys)
    rect(xleft=dims[1,1],ybottom=dims[2,1],
         xright=dims[1,2],ytop=dims[2,2],
         col= rgb(1,1,1, alpha=0.3))

    natural.cover <- crop(landcover.nat,  dims)
    only.natural.cover <-  natural.cover[!natural.cover@data$VegName
                                         %in% kinda.natural,]
    non.native.natural.cover <- natural.cover[natural.cover@data$VegName
                                              %in% kinda.natural,]

    colors <-
        add.alpha(rainbow((length(unique(natural.cover@data$VegName))
        +1 )),
                  0.5)
    names(colors) <- unique(natural.cover@data$VegName)

    for(veg.type in unique(natural.cover@data$VegName)){
        plot(natural.cover[natural.cover@data$VegName == veg.type,],
             col=colors[veg.type], add=TRUE)
    }

    points(new.hr, col="darkturquoise", pch=15)
    points(surveyed.hr, col="dodgerblue", pch=16)
    points(surveyed.control, col="navy", pch=17)

    points(new.hr, col="black", pch=0)
    points(surveyed.hr, col="black", pch=1)
    points(surveyed.control, col="black", pch=2)

    leg.names <- c("Riparian Forest Association",
                   "California Annual Grasslands Alliance",
                   "Upland Annual Grasslands & Forbs",
                   "Wetland Forbs Super Alliance",
                   "Valley Oak Alliance - Riparian",
                   "Flooded Deciduous Shrubland",
                   "Wet Meadow Grasses Super Alliance",
                   "Mixed Willow Super Alliance",
                   "Valley Oak Alliance", "Blackberry Super Alliance",
                   "Tamarisk Alliance", "Fremont Cottonwood/Willow",
                   "Fresh Water Marsh Alliance",
                   "Eucalyptus Alliance", "Blue Oak Alliance",
                   "Interior Live/Blue Oak", "Vernal Pool Complex")

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n")
    legend("left", col=colors, pch=16,
           legend=leg.names, ncol=2, bty="n")
    legend("left", col="black", pch=1,
           legend=leg.names, ncol=2, bty="n")

    only.2.cols <- add.alpha(c("red", "violetred"), 0.5)
    plot(sys)
    plot(non.native.natural.cover, col=only.2.cols[1], add=TRUE)
    plot(only.natural.cover, col=only.2.cols[2], add=TRUE)

    points(new.hr, col="darkturquoise", pch=15)
    points(surveyed.hr, col="dodgerblue", pch=16)
    points(surveyed.control, col="navy", pch=17)

    points(new.hr, col="black", pch=0)
    points(surveyed.hr, col="black", pch=1)
    points(surveyed.control, col="black", pch=2)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n")
    legend("left", col=only.2.cols, pch=16,
           legend=c("Predominantly non-native vegetation",
                    "Predominantly native vegetation"),
           ncol=1, bty="n")
    legend("left", col="black", pch=1,
           legend=c("Predominantly non-native vegetation",
                    "Predominantly native vegetation"),
           ncol=1, bty="n")


    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n")
    legend("left",
           col=c("darkturquoise","dodgerblue", "navy"),
           pch=c(15, 16, 17),
           legend=c("Unsurveyed hedgerow",
                    "Surveyed hedgerow",
                    "Surveyed field margin"),
           ncol=3, bty="n")
    legend("left", col="black",
           pch=c(0, 1, 2),
           legend=c("Unsurveyed hedgerow",
                    "Surveyed hedgerow",
                    "Surveyed field margin"),
           ncol=3, bty="n")

}


pdf.f(plotLandscapeCovers,
      file="../../../hedgerow_metacom_saved/map/covermap.pdf",
      height=10, width=8)


## looking at each type of vegetation on the landscape
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
