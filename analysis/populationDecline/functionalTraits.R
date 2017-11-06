rm(list=ls())
library(FD)
library(vegan)
setwd('~/Dropbox/hedgerow_metacom/analysis/populationDecline')
source("src/calcFuncUniqOrig.R")
load('../../../hedgerow_assembly/data/networks/specimens.Rdata')
traits <- read.csv("../../../hedgerow_assembly/data/traits.csv")

traits$mean.abun <- apply(traits, 1, function(x){
  mean(as.numeric(c(x["abun.net"], x["abun.pan"])),
       na.rm=TRUE)
})

bee.traits <- c("NestLoc", "Excavate", "Sociality","Lecty","MeanITD")
syr.traits <- c("LarvalHabitat", "LarvalDiet", "AdultDiet",
                "BodyLength")
shared.traits <- c("GenusSpecies", "degree", "d",
                   "temp.mean", "max.temp", "temp.range")
plant.traits <- c("GenusSpecies", "degree", "d")

bee.weights <- c(rep(1, 2), rep(1/3, 3), rep(1, length(bee.traits)))

bee.func <- calcFuncUniqOrig(traits,
                             traits.2.keep=c(bee.traits, shared.traits),
                             bee.nonbee="bee",
                             abund.col= "mean.abun.net",
                             weights=bee.weights)

syr.weights <- c(rep(1, 2), rep(1/3, 3), rep(1, length(syr.traits)))

syr.func <- calcFuncUniqOrig(traits,
                             traits.2.keep=c(syr.traits, shared.traits),
                             bee.nonbee="non-bee",
                             abund.col= "mean.abun.net",
                             weights=syr.weights)

plant.weights <- rep(1, 2)

plant.func <- calcFuncUniqOrig(traits,
                             traits.2.keep=plant.traits,
                             bee.nonbee="plant",
                             abund.col= "mean.abun.net",
                             weights=plant.weights)



func <- rbind(bee.func, syr.func, plant.func)

traits <- cbind(traits, func[match(traits$GenusSpecies,
                                 rownames(func)),])

write.csv(traits, file="../../data/traits.csv", row.names=FALSE)
