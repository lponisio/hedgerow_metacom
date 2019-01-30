## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')
source('plotting/src/predictIntervals.R')
source('plotting/src/CIplotting.R')
source('plotting/src/ZplotPanels.R')

source('src/initialize.R')
load(file=sprintf('saved/mods/zmets_drop_li_ht%s.Rdata', drop.li.ht))

ys <- c("degree", "betweenness", "closeness")
ylabs <- c("Degree", "Betweenness", "Closeness")

## ************************************************************
## pollinators
## ************************************************************
pol.specs <- list(year.sp, site.sp)

r.degree.dd <- expand.grid(r.degree=seq(from= min(year.sp$r.degree,
                                                  na.rm=TRUE),
                                        to= max(year.sp$r.degree,
                                                na.rm=TRUE),
                                        length=10),
                           MeanITD=mean(year.sp$MeanITD,
                                           na.rm=TRUE),
                           SiteStatus="all")

body.dd <- expand.grid(MeanITD=seq(from= min(year.sp$MeanITD,
                                                  na.rm=TRUE),
                                        to= max(year.sp$MeanITD,
                                                na.rm=TRUE),
                                        length=10),
                           r.degree=mean(year.sp$r.degree,
                                         na.rm=TRUE),
                       SiteStatus="all")

## nets.sites is within a site across years (time), nets.years is across
## sites within a year (space)
pp <- c("space", "time")
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

## ************************************************************
## sites
## ************************************************************

year.site <- year.site[year.site$Div !=
                                     min(year.site$Div),]

frd.dd <- expand.grid(Div=seq(from= min(year.site$Div,
                                        na.rm=TRUE),
                              to= max(year.site$Div,
                                      na.rm=TRUE),
                              length=10),
                      natArea=mean(year.site$natArea,
                                   na.rm=TRUE),
                      hrArea=mean(year.site$hrArea,
                                  na.rm=TRUE),
                      SiteStatus="all")

natArea.dd <- expand.grid(natArea=seq(from= min(year.site$natArea,
                                                na.rm=TRUE),
                                      to= max(year.site$natArea,
                                              na.rm=TRUE),
                                      length=10),
                          Div=mean(year.site$Div,
                                   na.rm=TRUE),
                          hrArea=mean(year.site$hrArea,
                                      na.rm=TRUE),
                          SiteStatus="all")


## nets.sites is within a site across years (time), nets.years is across
## sites within a year (space)
mods.site <- list(space=mod.years.site)
specs <- year.site

## plot pollinator and patch panels together
plot.panels.all()
