rm(list=ls())
## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/speciesLevel')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')

xvar <- c("r.degree", "BodyLength")
source('src/initialize.R')

## ************************************************************
## specialization
## ************************************************************

load(file=file.path(save.path, "mods/specmetrics.Rdata"))

ys <- c("proportional.generality.x", "degree.x", "k", "betweenness", "closeness")

ylabs <- c("Proportional Generality", "Degree", "K",
           "Betweenness centrality", "Closeness centrality")

r.degree.dd <- expand.grid(r.degree=seq(from= min(all.specs[[1]]$r.degree,
                                                  na.rm=TRUE),
                                        to= max(all.specs[[1]]$r.degree,
                                                na.rm=TRUE),
                                        length=10),
                           BodyLength=mean(all.specs[[1]]$BodyLength,
                                           na.rm=TRUE),
                           SiteStatus="all")

body.dd <- expand.grid(BodyLength=seq(from= min(all.specs[[1]]$BodyLength,
                                                  na.rm=TRUE),
                                        to= max(all.specs[[1]]$BodyLength,
                                                na.rm=TRUE),
                                        length=10),
                           r.degree=mean(all.specs[[1]]$r.degree,
                                         na.rm=TRUE),
                       SiteStatus="all")

## nets.sites is within a site across years (time), nets.years is across
## sites within a year (space)
pp <- c("space", "time")
mods <- list(mod.years, mod.sites)
names(mods) <- names(all.specs) <- pp


makePlots(pp=pp, xvar=xvar,
          ys=ys, dd=r.degree.dd,
          mods=mods, ylabs=ylabs,
          all.specs=all.specs,
          xs="r.degree",
          xlab= "Floral degree")

makePlots(pp=pp, xvar=xvar,
          ys=ys, dd=body.dd,
          mods=mods, ylabs=ylabs,
          all.specs=all.specs,
          xs="BodyLength",
          xlab= "Body size")


plot.panels()




