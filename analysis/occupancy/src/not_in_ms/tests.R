## setwd('~/Dropbox/hedgerow_metacom')
rm(list=ls())
setwd('analysis/occupancy')

## load data of interest
## random data checks

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




