rm(list=ls())
library(FD)
library(vegan)
setwd('~/Dropbox/collapse/analysis')
source("src/calcFuncUniqOrig.R")
load('../../hedgerow_assembly/data/networks/specimens.Rdata')
traits <- read.csv("~/Dropbox/hedgerow_assembly/data/traits.csv")

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

write.csv(traits, file="../data/traits.csv", row.names=FALSE)


## calcSpeciesCoords <- function(traitmat, abund, weights, ...){
##   ## calculating species' coordinates and the centroid
##   coords <- dbFD(traitmat, abund, w=weights,
##                  corr="cailliez", print.pco=T)$x.axes
##   ## now first calculating the centroid of the pollis present for each
##   ## site
##   centr <- list(NULL)
##   ## centr = averaged trait values for pollis present. 
##   for(i in 1:col(abund)){
##     pres <- which(abund[i,] > 0)
##     vec <- coords[as.logical(abund[i]),]
##     w <- abund[i, pres]
##     centr[[i]] <- apply(vec, 2, weighted.mean, w = w)
##   }
##   centr <- do.call(rbind, centr)
##   rownames(centr) <- rownames(abund)
##   return(list(coords=coords, centr=centr))
## }


## FD_measures <- function(coords, centr, abund){
##   ## adds a line for the centroid coodinates in each community
##   coords2  <- list(NULL)
##   for (i in 1:length(abund[,1])){
##     coords2[[i]] <- data.frame(rbind(coords, centr[i,]))
##     rownames(coords2[[i]])[dim(coords)[1]+1] <- "centr"
##   }
  
##   dists_centr   <- lapply(coords2, function(x){
##     d.poll<- as.matrix(dist(x, diag=TRUE, upper=TRUE))
##     for (i in 1:dim(d.poll)[1]) {d.poll[i,i] <- NA}
##     return(d.poll)
##   })
  
##   ## Originality: Distance to centroid of the species present
  
##   originality<- unlist(lapply(1:length(abund[,1]), function(i){
##     dists <- dists_centr[[i]]
##     pres <- which(abund[i,]>0)
##     dists[pres,dim(dists)[1]]
##   }))
  
##   ## Uniqueness: nearest neighbour among present species
  
##   uniqueness <- unlist(lapply(1:length(abund[,1]), function(i){
##     pres <- which(abund[i,]>0)
##     dists <- dists_centr[[i]][pres,pres]
##     uniq <- NULL
##     for (j in 1:length(pres)){uniq <- c(uniq, min(dists[j,], na.rm=T))}
##     return (uniq)
##   }))
  
##   ## also need to keep track of farm, species names and their abundances
  
##   measures <- do.call(rbind, lapply(1:length(abund[,1]), function(i){
##     pres <- which(abund[i,]>0)
##     abundance <- abund[i, pres]
##     sp <- names(abundance)
##     site <- rownames(abund)[i]
##     farm <-  rep(site, length(pres))
##     tab <- data.frame(farm=farm, sp=sp, abundance=as.numeric(abundance))
##     return (tab)
##   }))
  
##   ## merge all
##   measures$orig <- originality
##   measures$uniq <- uniqueness
  
##   return(as.data.frame(measures))
  
## }

## weight <- c(0.5, 0.5,1,0.5,0.5,1, rep(1/6, 6), rep(1,3))
