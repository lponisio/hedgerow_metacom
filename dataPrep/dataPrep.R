## setwd('~/Dropbox/hedgerow_metacom/')
rm(list=ls())
setwd('dataPrep')
source('src/misc.R')
source('src/prepNets.R')
source('src/specialization.R')
library(bipartite, quietly = TRUE)
library(fossil, quietly = TRUE)

## *******************************************************************
## load all the datasets
## *******************************************************************
load('~/Dropbox/hedgerow/data_sets/traditional/specimens-complete.RData')
load('~/Dropbox/hedgerow/data_sets/matrices/net/bee.syr.RData')
trait.dir <- '~/Dropbox/hedgerow/data_sets/traditional/functional_traits'
load('~/Dropbox/hedgerow/data_sets/misc/veg.Rdata')

## *******************************************************************
## subset data to only net samples and bees with IDs
## *******************************************************************

## minimum number of samples across year needed to be included in the
## analysis
sample.min <- 5
spec <- dd

## subset to net specimens
spec <- spec[spec$NetPan == 'net',]
## create species column
spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                                 spec$PlantSpecies))

## drop pollinators and plants without identifications
spec <-  spec[spec$PlantGenusSpecies != '',]

## subset to just bees
spec <- spec[spec$BeeNonbee == 'bee',]
spec <-  spec[spec$Species != '',]
spec$SiteStatus[spec$SiteStatus == 'restored'] <- 'maturing'
spec$Date <- as.Date(spec$Date)
spec$doy <- as.numeric(strftime(spec$Date, format='%j'))

occ.data <-
    try(load('~/Dropbox/hedgerow_metacom_saved/occupancy/all-5-0-2500-350.Rdata'),
        silent=TRUE)
if(!inherits(occ.data, "try-error")){
    keep.sites <- dimnames(model.input$data$X)[[1]]
    spec <- spec[spec$Site %in% keep.sites,]
}

## *******************************************************************
## create a giant network to calculate specialization
## *******************************************************************
agg.spec <- aggregate(list(abund=spec$GenusSpecies),
                      list(GenusSpecies=spec$GenusSpecies,
                           PlantGenusSpecies=spec$PlantGenusSpecies),
                      length)

nets.all <- samp2site.spp(agg.spec$PlantGenusSpecies,
                          agg.spec$GenusSpecies,
                          agg.spec$abund, FUN=sum)

all.specializations <- specieslevel(nets.all,
                       index=c('proportional generality',
                                            'degree',
                                            'd'))
## calculate rarified plant.pol degree
rare.plants.degree <- apply(nets.all, 1, chao1)
rare.pols.degree <- apply(nets.all, 2, chao1)

traits <- data.frame(GenusSpecies= unlist(sapply(all.specializations,
                                                 rownames)),
                     do.call(rbind, all.specializations))

traits$r.degree <-  rare.pols.degree[match(traits$GenusSpecies,
                                           names(rare.pols.degree))]
traits$r.degree[is.na(traits$r.degree)] <-
    rare.plants.degree[match(traits$GenusSpecies[is.na(traits$r.degree)],
                             names(rare.plants.degree))]

rownames(traits) <- NULL

## *******************************************************************
## drop forb and natural sites after network trait calculation
## *******************************************************************
samples <- apply(apply(mat, c(1,2), function(x) all(!is.na(x))), 1,
                 sum)
names(samples) <- gsub(":.+", "", names(samples))
spec <- spec[spec$Site %in% names(samples[samples > sample.min]),]

to.drop.status <- c('forb', 'natural')

## drop the sites outside of Yolo C because there is not natural data
sites.to.drop <- c('Martinez', 'PutahCreekForb', 'PutahCreek', 'RSlough')

spec <- spec[!spec$SiteStatus %in% to.drop.status,]
spec <- spec[!spec$Site %in% sites.to.drop,]


## print quantities of interest for manuscript
## total specimens
paste0("specimens=", nrow(spec))

## total species
paste0("species=",length(unique(spec$GenusSpecies)))

## sampling dates
paste0("samples=", length(unique(paste(spec$Site, spec$Date))))

## families and genera
paste0("family=",length(unique(spec$Family)))
paste0("genus=", length(unique(spec$Genus)))

## interactions
paste0("interactions=", length(unique(paste(spec$GenusSpecies,
                                            spec$PlantGenusSpecies))))


## *******************************************************************
## sampling table for manuscript
## *******************************************************************
site.table <- aggregate(list(Samples=spec$Date),
                        list(Year=spec$Year, Site=spec$Site),
                        function(x) length(unique(x)))

ms.table <- samp2site.spp(site=site.table$Site,
                          spp=site.table$Year,
                          abund=site.table$Samples,
                          FUN=sum)

write.csv(ms.table, file="../data/samples.csv")

BACI.site <- c('Barger', 'Butler', 'Hrdy', 'MullerB', 'Sperandio')
spec$SiteStatus[spec$Site %in% BACI.site] <- "maturing"

ms.table <- cbind(spec$SiteStatus[match(rownames(ms.table),
                                        spec$Site)], ms.table)

colnames(ms.table) <- c("Site type", colnames(ms.table)[-1])

ms.table <- ms.table[order(ms.table[, "Site type"], decreasing=TRUE),]

ms.table[, "Site type"][ms.table[, "Site type"] == "maturing" |
  ms.table[, "Site type"] == "mature" | ms.table[, "Site type"] == "restored"] <-
  "Hedgerow"

ms.table[, "Site type"][ms.table[, "Site type"] == "control"] <-
    "Field margin"

rownames(ms.table) <- paste0("HR", 1:nrow(ms.table))

rownames(ms.table)[ms.table[, "Site type"] == "Field margin"] <-
    paste0("FM", 1:sum(ms.table[, "Site type"] == "Field margin"))

write.table(ms.table,
            file="~/Dropbox/hedgerow_metacom_saved/occupancy/table/samples.txt",
            sep=" & ")

## *******************************************************************
## add various traits
## *******************************************************************
spec <- merge(spec, traits)
save(spec, file='../data/networks/allSpecimens.Rdata')
write.csv(traits, file='../data/traits.csv', row.names=FALSE)

## *******************************************************************
## create networks
## *******************************************************************

site.years <- aggregate(Year ~ Site, data=spec,
                         function(x) length(unique(x)))
sites.to.keep <- site.years$Site[site.years$Year >= 1]

spec.for.nets <- spec[spec$Site %in% sites.to.keep,]
nets <- breakNet(spec.dat=spec.for.nets, 'Site', 'Year')
## within a year, the network of shared species across sites
nets.year <- breakNetSpTemp(spec.dat=spec.for.nets, 'Year', 'Site')
## within a site, the network of shared species across years
nets.site <- breakNetSpTemp(spec.dat=spec.for.nets, 'Site', 'Year')

## save networks for each site, timeframe
f.path <- '../data/networks'
save(nets, file=file.path(f.path, 'all_networks_years.Rdata'))
save(nets.year, file=file.path(f.path, 'years_networks.Rdata'))
save(nets.site, file=file.path(f.path, 'sites_networks.Rdata'))

## *******************************************************************
## site-species level metric calculation
## *******************************************************************
specs <- calcSpec(nets, spec)
specs$closeness.log <- log(specs$closeness + 1)

## within a year, across sites
specs.years <- calcSpec(nets.year, spec)
## within a site across years
specs.site <- calcSpec(nets.site, spec)

specs.save.path <- '../analysis/networks/saved'
save(specs, specs.years, specs.site,
     file=file.path(specs.save.path, 'specs.Rdata'))


## drop super generalists
spec.for.nets <- spec.for.nets[!spec.for.nets$GenusSpecies %in%
                               c("Lasioglossum (Dialictus) incompletum",
                                 "Halictus tripartitus"), ]

nets <- breakNet(spec.dat=spec.for.nets, 'Site', 'Year')
## within a year, the network of shared species across sites
nets.year <- breakNetSpTemp(spec.dat=spec.for.nets, 'Year', 'Site')
## within a site, the network of shared species across years
nets.site <- breakNetSpTemp(spec.dat=spec.for.nets, 'Site', 'Year')

## save networks for each site, timeframe
f.path <- '../data/networks'
save(nets, file=file.path(f.path, 'all_networks_years_no_li_ht.Rdata'))
save(nets.year, file=file.path(f.path, 'years_networks_no_li_ht.Rdata'))
save(nets.site, file=file.path(f.path, 'sites_networks_no_li_ht.Rdata'))


## *******************************************************************
## veg data based on visitation
## *******************************************************************
abund.visits <- lapply(nets, function(x) apply(x, 1, sum))
plant.div.visits <- sapply(abund.visits, diversity)


div.visits <- lapply(nets, function(x) apply(x, 1, diversity))
plant.median.div.visits <- sapply(div.visits, median)

by.site$Site <- sapply(strsplit(by.site$Site, ":"), function(x) x[1])
by.site$div.visits <- plant.div.visits[match(paste(by.site$Site,
                                                   by.site$Year, sep="."),
                                             names(nets))]

by.site$median.div.visits <- plant.median.div.visits[match(paste(by.site$Site,
                                                   by.site$Year, sep="."),
                                             names(nets))]

save(by.site, file="../data/veg.Rdata")
