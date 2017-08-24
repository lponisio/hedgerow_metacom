rm(list=ls())
library(plyr)
library(lme4)
library(lmerTest)
setwd('~/Dropbox/collapse/analysis')
load('../../hedgerow_assembly/data/species/allsamples.Rdata')
load('../../hedgerow_assembly/data/networks/allSpecimens.Rdata')
traits <- read.csv("../data/traits.csv")
source("src/runModels.R")

drought.yrs <- 2013:2015

## pollinators 
pol <- calcDroughtvar(pol.mat, drought.yrs, traits, spec=spec)
plant <- calcDroughtvar(plant.mat, drought.yrs, traits, spec=spec)

ys <- c("score")
xvars <- c("d", "degree",
           "proportional.generality",
           "uniq",
           "originality")
fams <- rep("gaussian", length(xvars))

## full models
formulas <- list()
for(i in xvars){
  formulas[[i]] <- lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(i,
                                 "*drought", sep=""),
                           "(1|site)",
                           "(1|species)",
                           sep="+")))
  })
}
formulas <- unlist(formulas)

## pollinators, controls
dats <- list(pol, plant)
metrics <- as.character(unique(pol$metric))
site.type <- c("control", "mature")
combins <- length(dats)*length(metrics)*length(site.type)*length(xvars)
out.mods <- vector("list",
                   length=combins)
out.sums <- vector("list",
                   length=combins)
this.dats <- vector("list",
                    length=combins)

for (i in 1:length(dats)){
  for(j in metrics){
    print(j)
    for(k in site.type){
      print(k)
      this.dats[[i]][[j]][[k]] <- dats[[i]][!is.na(dats[[i]]$SiteStatus),]
      this.dats[[i]][[j]][[k]] <- this.dats[[i]][[j]][[k]][this.dats[[i]][[j]][[k]][,"SiteStatus"] == k,] 
      this.dats[[i]][[j]][[k]] <- this.dats[[i]][[j]][[k]][this.dats[[i]][[j]][[k]][,"metric"] == j,]
      out.mods[[i]][[j]][[k]] <- mapply(function(a, b)
                                   runMod(forms= a,
                                          fam= b,
                                          dats=this.dats[[i]][[j]][[k]]),
                                   a=formulas,
                                   b=fams,
                                   SIMPLIFY=FALSE)
      out.sums[[i]][[j]][[k]] <- lapply(out.mods[[i]][[j]][[k]], summary)
    }
  }
}


save(out.sums, out.mods, file="saved/mods.Rdata")

## pollinators
## [1] "mean.abund"
## [1] "control"
## 1) degree*drought, 2) prop.gen*drought 3) originality*drought
## [1] "mature"
## 1) degree*drought,  3) originality*drought

## [1] "occ"
## [1] "control"
## 1) degree*drought 2) prop.gen*drought, 3) uniq*drought
## [1] "mature"
## 1) d*drought, 2) degree*drought, 3) prop.gen*drought, 4)
## uniq*drought, 5) origin*drought

## [1] "cv.abund"
## [1] "control"
## 1) d*drought 2) degree*drought 3) originality*drought
## [1] "mature"
## no interactions

## [1] "cv.occ"
## [1] "control"
## 1) d*drought 2) degree*drought 3) originality*drought
## [1] "mature"
## 1) d*drought

## plants
## [1] "mean.abund"
## [1] "control"
## 1) degree*drought, 2) prop.gen*drought, 3) uniq*drought 
## [1] "mature"
## 1) degree*drought 2) prop.gen*drought 3) uniq*drought

## [1] "occ"
## [1] "control"
## 1) degree*drought 2) prop.gen*degree 3) uniq*drought 4)
## orig*drought
## [1] "mature"
## 1) prop.gen*drought 2) uniq*drought 3) orig*drought

## [1] "cv.abund"
## [1] "control"
## 1) degree*drought 2) uniq*drought
## [1] "mature"
## no interactions

## [1] "cv.occ"
## [1] "control"
## 1) degree*drought 2) uniq*drought
## [1] "mature"
## in interactions
