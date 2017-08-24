library(vegan)
library(dismo)
library(rgdal)
library(RColorBrewer)
library(maptools)
library(igraph)
source('src/misc.R')
source('plotting/src/makeNetworkFig.R')


load('../../data/networks/allSpecimens.Rdata')
load('../../data/networks/years_networks.Rdata')
load('../../data/networks/sites_networks.Rdata')
load('../../data/spatial/sitesAll.Rdata')
load('../../data/spatial/landcoverNat.Rdata')

spec$SiteStatus[spec$Site == 'Barger' |
                spec$Site == 'Butler'|
                spec$Site == 'Hrdy'|
                spec$Site == 'MullerB'|
                spec$Site == 'Sperandio'] <- 'mature'

spec$SiteStatus[spec$SiteStatus == "maturing"] <- "mature"


## prep map and site data
all.sites.pt <- spTransform(all.sites.pt,
                            CRS("+init=epsg:4326"))

landcover.nat <- spTransform(landcover.nat,
                            CRS("+init=epsg:4326"))

bbox.sites <- bbox(all.sites.pt)
bbox.nat <- bbox(landcover.nat)

bbox.all <- matrix(c(min(bbox.sites[1,1], bbox.nat[1,1]),
                   min(bbox.sites[2,1], bbox.nat[2,1]),
                   max(bbox.sites[1,2], bbox.nat[1,2]),
                   max(bbox.sites[2,2], bbox.nat[2,2])), ncol=2)
## + c(0, -0.001, -0.01, -0.15)

sys <- gmap(x=bbox.all,
            scale=2,
            type="satellite")
all.sites.pt <- spTransform(all.sites.pt,
                            sys@crs)
lat.long <- coordinates(all.sites.pt)
rownames(lat.long) <- all.sites.pt@data$df0


nets.year.sp <- lapply(nets.year, t)


landcover.nat <- spTransform(landcover.nat,
                            sys@crs)
