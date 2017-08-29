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

## bbox.all <- matrix(c(min(bbox.sites[1,1], bbox.nat[1,1]),
##                    min(bbox.sites[2,1], bbox.nat[2,1]),
##                    max(bbox.sites[1,2], bbox.nat[1,2]),
##                    max(bbox.sites[2,2], bbox.nat[2,2])), ncol=2)

bbox.all <- matrix(c(bbox.nat[1,1],
                     bbox.nat[2,1],
                     bbox.nat[1,2],
                     bbox.nat[2,2]), ncol=2)

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


load('../../../hedgerow_metacom_saved/occupancy/5-0-all.RData')
load('../../../hedgerow_metacom_saved/occupancy/runs/nimble_bees_noRain.Rdata')

nimble.sum <- ms.ms.nimble$model1$summary["nimble",,]

getSiteAve <- function(pattern, nimble.summary, model.input){
    nimble.summary <- nimble.summary[, grep(pattern,
                                            colnames(nimble.summary))]

    indexes <- gsub("gam.site.mean\\[|phi.site.mean\\[","", colnames(nimble.summary))
    site.num <- sapply(strsplit(indexes, ","), function(x) x[1])

    nimble.site.ave <- tapply(nimble.summary["mean",], site.num, mean)

    sites <- dimnames(model.input$data$X)[[1]]

    names(nimble.site.ave) <- sites
    return(nimble.site.ave)
}

phi.site.ave <- getSiteAve("phi.site.mean", nimble.sum, model.input)
gam.site.ave <- getSiteAve("gam.site.mean", nimble.sum, model.input)


phi.gam <- list(phi.site.ave, gam.site.ave*2)
