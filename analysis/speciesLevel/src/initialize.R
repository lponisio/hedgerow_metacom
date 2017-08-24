library(igraph)
library(bipartite)
library(lme4)
library(lmerTest)
library(RColorBrewer)
source('src/misc.R')
source("src/specialization.R")
source("src/prepNets.R")

traits <- read.csv("../../data/traits.csv")
load('../../../hedgerow_metacom/data/networks/allSpecimens.Rdata')
load('../../../hedgerow_metacom/data/networks/all_networks_years.Rdata')


save.path <- 'saved'

load(file='saved/specs.Rdata')

## trait data
traits <- read.csv("../../data/traits.csv")
traits$Nestuse <- paste(traits$NestLoc, traits$Excavate)
traits$Nestuse[traits$Nestuse == "NA NA"] <- NA
traits$Nestuse[traits$Nestuse == "NA rent"] <- NA
traits$BodyLength[!is.na(traits$MeanITD)] <-
    traits$MeanITD[!is.na(traits$MeanITD)]

specs.site <- merge(specs.site, traits, by="GenusSpecies")
specs.years <- merge(specs.years, traits, by="GenusSpecies")

specs.years.pol <- specs.years[specs.years$speciesType ==
                               "pollinator",]
specs.site.pol <- specs.site[specs.site$speciesType ==
                             "pollinator",]

specs.site.pol[, xvar] <- apply(specs.site.pol[, xvar], 2,
                                 standardize)
specs.years.pol[, xvar] <- apply(specs.years.pol[, xvar], 2,
                                 standardize)

all.specs <- list(specs.years.pol, specs.site.pol)


