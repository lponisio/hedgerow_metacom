library(igraph, quietly = TRUE)
library(bipartite,  quietly = TRUE)
library(lme4,  quietly = TRUE)
library(lmerTest,  quietly = TRUE)
library(RColorBrewer,  quietly = TRUE)
source('src/misc.R')
source("src/specialization.R")
source("src/prepNets.R")

fig.path <- "../../../hedgerow_metacom_saved/networks/figures"
checkDirExists(fig.path)

load('../../../hedgerow_metacom/data/networks/allSpecimens.Rdata')
load('../../../hedgerow_metacom/data/networks/all_networks_years.Rdata')


if(!exists("passed.args")){
    passed.args <- commandArgs(trailingOnly=TRUE)
}
print(passed.args)
drop.li.ht <- as.logical(passed.args[1])
natural.decay <- passed.args[2]
HR.decay <- passed.args[3]

xvar.species <- c("r.degree", "MeanITD")
xvar.site <- c("Div", "natArea", "hrArea")


## species traits
traits <- read.csv("../../data/traits.csv")

## hr proximity
load('../../data/spatial/hrcover_decay.Rdata')
## remnant proximity
load('../../data/spatial/natcover_decay_yolo.Rdata')
## veg data
load('../../data/veg.Rdata')
by.site$Site <- gsub(":.*", "", by.site$Site)
by.site$SiteYear <- paste(by.site$Site, by.site$Year)
by.site$natArea <- nat.area.sum[, natural.decay][match(by.site$Site,
                                                       rownames(nat.area.sum))]

by.site$hrArea <- hr.area.sum[, HR.decay][match(by.site$Site,
                                                rownames(nat.area.sum))]

save.path <- 'saved'
if(!drop.li.ht){
    load(file=file.path(save.path, 'specs.Rdata'))
}else if(drop.li.ht){
    load(file=file.path(save.path, 'specs_no_li_ht.Rdata'))
}

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
traits$MeanITD[!is.na(traits$MeanITD)] <-
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


save(specs.years.pol, specs.site.pol, specs.site.site,
     specs.years.site,
     file=file.path(save.path, 'specs_site_pol.Rdata'))
