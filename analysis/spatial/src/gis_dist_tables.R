rm(list=ls())
setwd('~/Dropbox/baci')
source('data/make/src/gis.R')
library(parallel)
library(maptools)
library(rgdal)
library(rgeos)
library(spatstat)
## ************************************************************

## ************************************************************
## load and prepare data
projection <- CRS('+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0')
## load landcover layers
load.and.merge <- function(year) {
  load(sprintf('../hedgerow/data_sets/gis/landcover/digitize%d.RData', year))
  dd <- spTransform(dd, projection)
  unionSpatialPolygons(SpP=dd, IDs=dd@data$LandCover)
}
dd.2006 <- load.and.merge(2006)
dd.2012 <- load.and.merge(2012)

## load sites and subset down to just BACI
load('../hedgerow/data_sets/gis/hedgerows.RData')
dd.h <- spTransform(dd.h, projection)

baci <- c('Barger',
          'Butler',
          'Hrdy',
          'MullerB',
          'Sperandio',
          'BC2',
          'Chamberlain',
          'DQU',
          'Gregory',
          'H16',
          'HC1',
          'MC1',
          'Spiva',
          'Turkovich',
          'USS',
          'BC1',
          'Gnoss',
          'Rock',
          'GC1')

dd.h <- dd.h[dd.h@data$site %in% baci,]

radii <- round(exp(seq(from=log(10), to=log(1500), length=20)))

dist.tables.2006 <- make.distance.table(dd.lc=dd.2006,
                                        year=2006,
                                        radii=radii,
                                        num.cores=8)

dist.tables.2012 <- make.distance.table(dd.lc=dd.2012,
                                        year=2012,
                                        radii=radii,
                                        num.cores=8)

dist.tables <- list('2006'=dist.tables.2006,
                    '2012'=dist.tables.2012)

save(dist.tables,
     file='~/Desktop/dist.tables.RData')


save(dist.tables,
     file='data/saved/dist.tables.RData')
## ************************************************************



## ## identify sites with cover data for 2006
## plot(dd.2006)
## plot(dd.h, col='red', add=TRUE)

## ## identify sites with cover data for 2012
## plot(dd.2012)
## plot(dd.h, col='red', add=TRUE)
## plot(dd.h[dd.h@data$site %in% c('BC1', 'BC2', 'GC1', 'HC1',
##                                 'Gnoss', 'Rock'),],
##      col='yellow', add=TRUE)
