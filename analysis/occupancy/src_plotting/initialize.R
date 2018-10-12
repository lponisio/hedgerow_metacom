library(vegan)
library(dismo)
library(rgdal)
library(RColorBrewer)
library(maptools)
## library(maps)
library(igraph)
library(bipartite)
library(rgeos)
source('src/misc.R')

save.dir <- "../../../hedgerow_metacom_saved/occupancy"
checkDirExists(save.dir)
checkDirExists(file.path(save.dir, "figures/networks"))

load('../../data/networks/allSpecimens.Rdata')
load('../../data/networks/years_networks.Rdata')
load('../../data/networks/sites_networks.Rdata')
load('../../data/spatial/sitesAll.Rdata')
load('../../data/spatial/landcoverNatPatches.Rdata')

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

bbox.all <- matrix(c(bbox.sites[1,1],
                     bbox.sites[2,1],
                     bbox.sites[1,2],
                     bbox.sites[2,2]),
                   ncol=2)

makeSys <- function(){
    sys <- gmap(x=bbox.all,
                    scale=2,
                    type="satellite", zoom=11)
}

## use the phi and gamma data from the occupancy model to make the
## figures
## load('../../../hedgerow_metacom_saved/occupancy/runs/all_nimble_bees_2500_350.Rdata')

## getSiteAve <- function(pattern, nimble.summary, model.input){
##     nimble.summary <- nimble.summary[, grep(pattern,
##                                             colnames(nimble.summary))]

##     indexes <- gsub("gam.site.mean\\[|phi.site.mean\\[","", colnames(nimble.summary))
##     site.num <- sapply(strsplit(indexes, ","), function(x) x[1])

##     nimble.site.ave <- tapply(nimble.summary["mean",], site.num, mean)

##     sites <- dimnames(model.input$data$X)[[1]]

##     names(nimble.site.ave) <- sites
##     return(nimble.site.ave)
## }

## phi.site.ave <- getSiteAve("phi.site.mean", nimble.sum, model.input)
## gam.site.ave <- getSiteAve("gam.site.mean", nimble.sum, model.input)


## phi.gam <- list(phi.site.ave, gam.site.ave*4)


## use the centrality scores

## maybe not very interesting given was not significant?
## load('../speciesLevel/saved/specs_site_pol.Rdata')


## site.between.mean <- tapply(specs.years.site$betweenness,
##                             specs.years.site$GenusSpecies, mean)*20




