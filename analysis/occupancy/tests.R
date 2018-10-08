## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')
args <- commandArgs(trailingOnly=TRUE)
args <- c("allInt","2500", "350", "filtering","all",1e2)
source('src/initialize.R')

## ************************************************************
## prep data
## ************************************************************

model.input <- prepOccModelInput(nzero=0,
                    threshold=2,
                    save.dir=save.dir,
                    spec,
                    sr.sched,
                    all.traits,
                    col.name.trait1 = "r.degree",
                    col.name.trait2 = "BodyLength",
                    HRarea= hr.area.sum, #sum.dist.area, ##spstats
                    natural.mat=nat.area.sum, ## natural
                    natural.decay=natural.decay,
                    HR.decay=HR.decay,
                    veg=by.site, #raw.flower.data,
                    load.inits=FALSE,
                    model.type=include.int,
                    col.name.div.type = "Div",## div.visits, Div
                    raw.flower=FALSE,
                    drop.li.ht=FALSE,
                    only.li.ht=FALSE)



## ************************************************************
## colinearily?
## ************************************************************

cor.test(model.input$data$natural, model.input$data$HRarea)

apply(model.input$data$fra, 2,
      function(x) cor.test(x, model.input$data$natural))

apply(model.input$data$fra, 2,
      function(x) cor.test(x, model.input$data$HRarea))

## nothing raises any red flags


## ************************************************************
## correlation between visit based flower diversity and pollinator
## abund and div?
## ************************************************************
library(vegan)
cor.test(by.site$Div, by.site$div.visits)

## calculate site-level statistics
calc.site.level <- function(dats, abund.type="sum"){
  rich <- length(unique(dats))
  dats <- as.character(dats)
  div <-  diversity(table(dats), index="shannon")
  if(abund.type == "median"){
    abund <- median(table(dats))
  } else {
    abund <- length(dats)
  }
  return(c(rich=rich, abund=abund, div=div))
}


by.site.pols <- aggregate(spec$GenusSpecies,
                       list(Site=spec$Site,
                            Year=spec$Year),
                       function(x) calc.site.level(x))

by.site <- merge(by.site, by.site.pols, by=c("Site", "Year"))

cor.test(by.site$x[,"rich"], by.site$div.visits)
cor.test(by.site$x[,"div"], by.site$div.visits)
cor.test(by.site$x[,"abund"], by.site$div.visits)

cor.test(by.site$x[,"rich"], by.site$Div)
cor.test(by.site$x[,"div"], by.site$Div)
cor.test(by.site$x[,"abund"], by.site$Div)




