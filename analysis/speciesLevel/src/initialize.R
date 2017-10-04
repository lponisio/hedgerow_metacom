library(igraph)
library(bipartite)
library(lme4)
library(lmerTest)
library(RColorBrewer)
source('src/misc.R')
source("src/specialization.R")
source("src/prepNets.R")

load('../../../hedgerow_metacom/data/networks/allSpecimens.Rdata')
load('../../../hedgerow_metacom/data/networks/all_networks_years.Rdata')

## species traits
traits <- read.csv("../../data/traits.csv")

## area weighted by log distance
load('../../data/spatial/HRareaDist.Rdata')
## area of different landcovers
load('../../data/spatial/natcover_decay_yolo.Rdata')
## veg data
load('../../data/veg.Rdata')
by.site$Site <- gsub(":.*", "", by.site$Site)
by.site$SiteYear <- paste(by.site$Site, by.site$Year)
by.site$natArea <- nat.area.sum[, natural.decay][match(by.site$Site,
                                                       rownames(nat.area.sum))]

by.site$hrArea <- sum.dist.area[match(by.site$Site,
                                                       names(sum.dist.area))]

save.path <- 'saved'
load(file='saved/specs.Rdata')

## specs.year is within a year across sites
## specs.site is within a site across years

## by site
specs.years.site <- specs.years[specs.years$speciesType ==
                               "plant",]
specs.site.site <- specs.site[specs.site$speciesType ==
                             "plant",]

specs.years.site$SiteYear <- paste(specs.years.site$GenusSpecies,
                                   specs.years.site$Site)

specs.site.site$SiteYear <- paste(specs.site.site$Site,
                                   specs.site.site$GenusSpecies)

specs.site.site <- merge(specs.site.site, by.site, by="SiteYear")
specs.years.site <- merge(specs.years.site, by.site, by="SiteYear")

specs.site.site[, xvar.site] <- apply(specs.site.site[, xvar.site], 2,
                                 standardize)
specs.years.site[, xvar.site] <- apply(specs.years.site[, xvar.site], 2,
                                       standardize)

specs.years.site$Site <- specs.years.site$Site.x
specs.site.site$Site <- specs.site.site$Site.x

## use the name naming as the pollinator dataset
specs.years.site$proportional.generality.x <-
    specs.years.site$proportional.generality
specs.site.site$proportional.generality.x <-
    specs.site.site$proportional.generality

specs.years.site$degree.x <- specs.years.site$degree
specs.site.site$degree.x <- specs.site.site$degree

## trait data
traits$Nestuse <- paste(traits$NestLoc, traits$Excavate)
traits$Nestuse[traits$Nestuse == "NA NA"] <- NA
traits$Nestuse[traits$Nestuse == "NA rent"] <- NA
traits$BodyLength[!is.na(traits$MeanITD)] <-
    traits$MeanITD[!is.na(traits$MeanITD)]

specs.years.pol <- specs.years[specs.years$speciesType ==
                               "pollinator",]
specs.site.pol <- specs.site[specs.site$speciesType ==
                             "pollinator",]

specs.site.pol <- merge(specs.site.pol, traits, by="GenusSpecies")
specs.years.pol <- merge(specs.years.pol, traits, by="GenusSpecies")

specs.site.pol[, xvar.species] <- apply(specs.site.pol[, xvar.species], 2,
                                 standardize)
specs.years.pol[, xvar.species] <- apply(specs.years.pol[, xvar.species], 2,
                                 standardize)


