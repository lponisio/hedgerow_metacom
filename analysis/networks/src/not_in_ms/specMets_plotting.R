rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/speciesLevel')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')
source('src/initialize.R')

## ************************************************************
## specialization
## ************************************************************

load(file=file.path(save.path, 'specs.Rdata'))
load(file=file.path(save.path, "mods/specs_ypr.Rdata"))

ylabs <- c("Proportional Generality", "Specialization (d')", "Degree",
           "Betweenness", "Closeness" )

dd <- expand.grid(ypr=seq(from= min(specs$ypr, na.rm=TRUE),
                          to= max(specs$ypr, na.rm=TRUE),
                          length=10))

pp <- c("plants", "pols")
mods <- list(mod.pols, mod.plants)
names(mods) <- pp

for(j in pp){
  for(i in 1:length(ys)){
    dd.ypr <- cbind(dd, 0)
    colnames(dd.ypr) <- c("ypr", ys[i])

      ypr.pi <- predict.int(mod= mods[[j]][[i]],
                            dd=dd.ypr,
                            y=ys[i],
                            family="gaussian")

    plot.predict.ypr(new.dd=ypr.pi,
                     ylabel=ylabs[i],
                     dats=specs,
                     y1=ys[i],
                     extinction.method=j,
                     agg.col="GenusSpecies")
  }
}


## closeness only

ypr.pi.pol <- predict.int(mod= mods$pol$closeness.log,
                      dd=dd.ypr,
                      y="closeness.log",
                      family="gaussian")


ypr.pi.plant <- predict.int(mod= mods$plant$closeness.log,
                      dd=dd.ypr,
                      y="closeness.log",
                      family="gaussian")

plot.panels()
