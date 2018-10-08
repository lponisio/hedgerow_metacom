library(sp)
library(maptools)
library(rgdal)
library(raster)
library(rgeos)

source('src/CalcSpStats.R')
source('src/getNatural.R')
source('src/makeDistanceTable.R')
load('../../data/spatial/hr.Rdata')
load('../../data/spatial/sitesAll.Rdata')
load('../../data/spatial/sitesAllLines.Rdata')

## area of different landcovers from kremen lab
load('../../data/spatial/notcreated/distTables.Rdata')
## area of landcover from yolo county
load('../../data/spatial/landcoverNat.Rdata')
## specimen data
load('../../data/networks/allSpecimens.Rdata')
