## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/spatial')
library(sp)
library(maptools)

source('src/CalcSpStats.R')
source('src/getNatural.R')
source('src/makeDistanceTable.R')
load('../../data/spatial/hr.Rdata')
load('../../data/spatial/sitesAll.Rdata')
load('../../data/spatial/sitesAllLines.Rdata')
## area of different landcovers from kremen lab
load('../../data/spatial/distTables.Rdata')
## area of landcover from yolo county
load('../../data/spatial/landcoverNat.Rdata')
## specimen data
load('../../data/networks/allSpecimens.Rdata')

## *************************************************************
## hedgerow area within buffers
## *************************************************************
buff <- seq(500, 2000, by=500)
rast.hr <- raster(hr)
nsite <- nrow(all.sites.pt@data)

spstats <- t(sapply(1:nsite, calcSpStats, d = buff,
                  plt = all.sites.pt, plt.name = 'SITE',
                  rast = hr))

colnames(spstats) <- paste("d", buff, sep="")
rownames(spstats) <- all.sites.pt@data$df0

save(spstats, file="../../data/spatial/HRarea.Rdata")

## *************************************************************
## natural using kremen lab digitized data, decay method
## *************************************************************
decays <- c(350, 1000, 2500)

natural.land <- getNatural(spec, dist.tables, decays)
save(natural.land,
     file="../../data/spatial/natcover.Rdata")

samp.site.nat.krem <- natural.land[natural.land$Site %in%
                                   spec$Site,]

## *************************************************************
## hedgerow area using decay method
## *************************************************************
radii <- round(exp(seq(from=log(10), to=log(1500), length=20)))
all.sites.lines@data$type <- "hedgerow"

hr.area <- makeDistanceTable(dd.lc=all.sites.lines,
                              dd.h=all.sites.pt,
                                        radii=radii,
                              num.cores=1,
                             type.col="type",
                             area.function=gLength,
                             drop.self=TRUE)
hr.area <- hr.area[!sapply(hr.area, is.null)]

hr.area.sum <- sapply(decays, function(x) {
    sapply(hr.area, applyDecay, decay=x)
})
colnames(hr.area.sum) <- decays

save(hr.area.sum,
     file="../../data/spatial/hrcover_decay.Rdata")

## *************************************************************
## natural using yolo county data source, decay method
## *************************************************************

nat.area <- makeDistanceTable(dd.lc=landcover.nat,
                              dd.h=all.sites.pt,
                                        radii=radii,
                              num.cores=2,
                              type.col="type")

nat.area <- lapply(nat.area, function(x){
    if(inherits(x, "try-error")){
        x <- NULL
    } else{
        x
    }
})

nat.area <- nat.area[!sapply(nat.area, is.null)]
save(nat.area, file="../../data/spatial/natcover_yolo.Rdata")

nat.area.sum <- sapply(decays, function(x) {
    sapply(nat.area, applyDecay, decay=x)
})
colnames(nat.area.sum) <- decays

only.nat.sites <- crop(all.sites.pt, extent(landcover.nat))
in.nat <- all.sites.pt@data$df0[all.sites.pt@data$df0 %in%
                                only.nat.sites@data$df0]
nat.area.sum[!rownames(nat.area.sum) %in% in.nat] <- NA
hist(nat.area.sum)
save(nat.area.sum,
     file="../../data/spatial/natcover_decay_yolo.Rdata")

samp.site.nat <- nat.area.sum[rownames(nat.area.sum) %in% spec$Site,]

## compare kremen and yolo data layers

samp.site.nat <- as.data.frame(samp.site.nat)
samp.site.nat$Site <- rownames(samp.site.nat)

all.decays <- merge(samp.site.nat.krem, samp.site.nat, by="Site")

plot(all.decays$natural350[all.decays$Year == 2006],
     all.decays$'350'[all.decays$Year == 2006],
     xlab= "Kremen",
     ylab="Yolo", pch=16)

text(all.decays$natural350[all.decays$Year == 2006],
     all.decays$'350'[all.decays$Year == 2006],
     labels=all.decays$Site[all.decays$Year == 2006],
     cex=0.7, pos=3)
abline(a=0, b=1)
