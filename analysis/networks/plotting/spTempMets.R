## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/plotPanels.R')

source('src/initialize.R')
load(file=sprintf('saved/mods/specmetrics_drop_li_ht%s.Rdata', drop.li.ht))

ys <- c("k", "weighted.betweenness")

ylabs <- c("Degree (z)", "Weighted centrality")

## ************************************************************
## pollinators
## ************************************************************

pol.specs <- list(specs.site.pol, specs.years.pol)

r.degree.dd <- expand.grid(r.degree=seq(from= min(pol.specs[[1]]$r.degree,
                                                  na.rm=TRUE),
                                        to= max(pol.specs[[1]]$r.degree,
                                                na.rm=TRUE),
                                        length=10),
                           MeanITD=mean(pol.specs[[1]]$MeanITD,
                                           na.rm=TRUE),
                           SiteStatus="all")

body.dd <- expand.grid(MeanITD=seq(from= min(pol.specs[[1]]$MeanITD,
                                                  na.rm=TRUE),
                                        to= max(pol.specs[[1]]$MeanITD,
                                                na.rm=TRUE),
                                        length=10),
                           r.degree=mean(pol.specs[[1]]$r.degree,
                                         na.rm=TRUE),
                       SiteStatus="all")


## nets.sites is within a site across years (time), nets.years is across
## sites within a year (space)
pp <- c("space", "time")
pol.specs <- list(specs.site.pol, specs.years.pol)
pol.mods <- list(mod.years.pol, mod.sites.pol)
names(pol.mods) <- names(pol.specs) <- pp

makePlots(pp=pp, xvar=xvar.species,
          ys=ys, dd=r.degree.dd,
          mods=pol.mods, ylabs=ylabs,
          all.specs=pol.specs,
          xs="r.degree",
          xlab= "Floral diet breadth")

makePlots(pp=pp, xvar=xvar.species,
          ys=ys, dd=body.dd,
          mods=pol.mods, ylabs=ylabs,
          all.specs=pol.specs,
          xs="MeanITD",
          xlab= "Body size")

plot.panels()


## ************************************************************
## sites
## ************************************************************

## drop extrema to make plotting nicer
specs.years.site <- specs.years.site[specs.years.site$Div !=
                                     min(specs.years.site$Div),]

frd.dd <- expand.grid(Div=seq(from= min(specs.years.site$Div,
                                        na.rm=TRUE),
                              to= max(specs.years.site$Div,
                                      na.rm=TRUE),
                              length=10),
                      natArea=mean(specs.years.site$natArea,
                                   na.rm=TRUE),
                      hrArea=mean(specs.years.site$hrArea,
                                  na.rm=TRUE),
                      SiteStatus="all")

natArea.dd <- expand.grid(natArea=seq(from= min(specs.years.site$natArea,
                                                na.rm=TRUE),
                                      to= max(specs.years.site$natArea,
                                              na.rm=TRUE),
                                      length=10),
                          Div=mean(specs.years.site$Div,
                                   na.rm=TRUE),
                          hrArea=mean(specs.years.site$hrArea,
                                      na.rm=TRUE),
                          SiteStatus="all")


## nets.sites is within a site across years (time), nets.years is across
## sites within a year (space)
mods.site <- list(space=mod.years.site)
specs <- specs.years.site
plot.panels.sites()


plot.panels.all()
