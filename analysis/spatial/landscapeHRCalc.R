rm(list=ls())
## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/spatial')
library(sp)
library(maptools)

source('src/CalcSpStats.R')
source('src/getNatural.R')
source('src/makeDistanceTable.R')
load('../../data/spatial/hr.Rdata')
load('../../data/spatial/sitesAll.Rdata')
## area of different landcovers from kremen lab
load('../../data/spatial/distTables.Rdata')
## area of landcover from yolo county
load('../../data/spatial/landcoverNat.Rdata')
load('../../data/spatial/landcoverNat_raster.Rdata')
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


## taking a look at the sampled sites natural habitat area
samp.site <- spstats[rownames(spstats) %in% spec$Site,]
write.csv(samp.site, file="../../data/spatial/hrarea.csv")

## *************************************************************
## natural using yolo county data source, buffers
## *************************************************************

## spstats.landcover <- t(sapply(1:nsite, calcSpStats, d = buff,
##                   plt = all.sites.pt, plt.name = 'SITE',
##                   rast = landcover.r))


## colnames(spstats) <- paste("d", buff, sep="")
## rownames(spstats) <- all.sites.pt@data$df0

## *************************************************************
## natural using kremen lab digitized data, decay method
## *************************************************************
decays <- c(350, 1000, 2500)

natural.land <- getNatural(spec, dist.tables, decays)
save(natural.land,
     file="../../data/spatial/natcover.Rdata")

samp.site.nat.krem <- natural.land[natural.land$Site %in%
                                   spec$Site,]

write.csv(samp.site.nat.krem, file="../../data/spatial/natarea_kremen.csv")


## *************************************************************
## natural using yolo county data source, decay method
## *************************************************************
landcover.nat@data$type <- "natural"

radii <- round(exp(seq(from=log(10), to=log(1500), length=20)))

nat.area <- makeDistanceTable(dd.lc=landcover.nat,
                              dd.h=all.sites.pt,
                                        radii=radii,
                              num.cores=1,
                              type.col="type")

save(nat.area, file="../../data/spatial/natcover_yolo.Rdata")

applyDecay <- function(dist.tab, decay){
    dist.tab <- dist.tab[[1]]
    areas <- c(dist.tab$area[1], diff(dist.tab$area))
    rads <- c(0, dist.tab$radius[1:(nrow(dist.tab)-1)])
    decay.area <- sum(exp(-(rads-0)^2/(2*decay^2)) * areas)
    return(decay.area)
}

nat.area.sum <- sapply(decays, function(x) {
    sapply(nat.area, applyDecay, decay=x)
})
colnames(nat.area.sum) <- decays

only.nat.sites <- crop(all.sites.pt, extent(landcover.nat))
in.nat <- all.sites.pt@data$df0[all.sites.pt@data$df0 %in%
                                only.nat.sites@data$df0]
nat.area.sum[!rownames(nat.area.sum) %in% in.nat] <- NA


save(nat.area.sum, file="../../data/spatial/natcover_decay_yolo.Rdata")



samp.site.nat <- nat.area.sum[rownames(nat.area.sum) %in% spec$Site,]
write.csv(samp.site, file="../../data/spatial/hrarea.csv")


write.csv(samp.site.nat,
          file="../../data/spatial/natarea_yolo.csv")

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
