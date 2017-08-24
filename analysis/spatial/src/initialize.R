library(sp)
library(raster)
library(spatstat)
library(maptools)
library(SDMTools)
library(rgdal)
source('src/CalcSpStats.R')

data.path <- '../../../data/spatial/'
save.path <- '../../../saved/spatial/'

## Read in data
farm <- na.omit(read.csv(file.path(data.path, 'landscape.csv')))
tanz.farm <- farm[farm$country == "TZA",]
gha.farm <- farm[farm$country == "GHA",]
ug.farm <- farm[farm$country == "UGA",]

tanz <- readGDAL(file.path(data.path, 'tanz/tz_tm1-57palglobdem_landcover.dat'))
tanz <- raster(tanz)

gha <- readGDAL(file.path(data.path, 'ghana/gh_tm1-57palglobdem_landcover.dat'))
gha <- raster(gha)

ug <- readGDAL(file.path(data.path, 'uganda/ug_tm1-57palglobdem_landcover.dat'))
ug <- raster(ug)
