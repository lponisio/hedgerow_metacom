rm(list=ls())
## setwd('~/Dropbox/hedgerow_metacom/')
setwd('dataPrep')
source('src/misc.R')
source('src/prepNets.R')
source('src/specialization.R')
library(bipartite)
library(fossil)

## load all the datasets
load('~/Dropbox/hedgerow/data_sets/traditional/specimens-complete.RData')
load('~/Dropbox/hedgerow/data_sets/matrices/net/bee.syr.RData')
trait.dir <- '~/Dropbox/hedgerow/data_sets/traditional/functional_traits'
bee.trait <-
    read.csv(file.path(trait.dir, 'bee.csv'),
             row.names=1)
both.trait <-
    read.csv(file.path(trait.dir, 'bee.syr.csv'))

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
    try(load('~/Dropbox/hedgerow_metacom_saved/occupancy/5-0-all.RData'),
        silent=TRUE)
if(!inherits(occ.data, "try-error")){
    keep.sites <- dimnames(model.input$data$X)[[1]]
    spec <- spec[spec$Site %in% keep.sites,]
}


## *************************************************
## create a giant network to calculate specialization
## *************************************************
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

## *************************************************
## drop forb and natural sites
## *************************************************
samples <- apply(apply(mat, c(1,2), function(x) all(!is.na(x))), 1,
                 sum)
names(samples) <- gsub(":.+", "", names(samples))

spec <- spec[spec$Site %in% names(samples[samples > sample.min]),]


to.drop.status <- c('forb', 'natural')

## drop the sites outside of Yolo C because there is not natural data
sites.to.drop <- c('Martinez', 'PutahCreekForb', 'PutahCreek', 'RSlough')

spec <- spec[!spec$SiteStatus %in% to.drop.status,]
spec <- spec[!spec$Site %in% sites.to.drop,]

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
paste0("interactions=", length(unique(paste(spec$GenusSpecies, spec$PlantGenusSpecies))))


## *************************************************
## sampling table for manuscript
## *************************************************

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

write.table(ms.table, file="~/Dropbox/hedgerow_metacom_saved/occupancy/table/samples.txt",
            sep=" & ")



## *************************************************
## thermal traits
## *************************************************
spec$AverageTemp <- apply(spec, 1, function(x){
    mean(as.numeric(c(x['TempStart'], x['TempEnd'])),
         na.rm=TRUE)
})

temp.tol <- do.call(rbind, tapply(spec$AverageTemp, spec$GenusSpecies,
                                  function(x){
                                      temp.mean <- mean(x, na.rm=TRUE)
                                      temp.range <- range(x, na.rm=TRUE)
                                      max.temp <- temp.range[2]
                                      temp.range <- temp.range[2] - temp.range[1]
                                      return(c(temp.mean=temp.mean,
                                               max.temp=max.temp,
                                               temp.range=temp.range))
                                  }))
temp.tol <- as.data.frame(temp.tol)
temp.tol$GenusSpecies <- rownames(temp.tol)
rownames(temp.tol) <- NULL

traits <- merge(traits, temp.tol, all.x=TRUE)


## *************************************************
## phenology traits
## *************************************************

doy.tol <- do.call(rbind, tapply(spec$doy, spec$GenusSpecies,
                                 function(x){
                                     doy.mean <- mean(x, na.rm=TRUE)
                                     doy.range <- range(x, na.rm=TRUE)
                                     max.doy <- doy.range[2]
                                     doy.range <- (doy.range[2] -
                                                   doy.range[1]) +1
                                     return(c(doy.mean=doy.mean,
                                              max.doy=max.doy,
                                              doy.range=doy.range))
                                 }))
doy.tol <- as.data.frame(doy.tol)
doy.tol$GenusSpecies <- rownames(doy.tol)
rownames(doy.tol) <- NULL

traits <- merge(traits, doy.tol, all.x=TRUE)


## *************************************************
## add various traits
## *************************************************
## specialization
spec$d <- traits$d[match(spec$GenusSpecies, traits$GenusSpecies)]
spec$degree <- traits$degree[match(spec$GenusSpecies,
                                   traits$GenusSpecies)]
spec$plant.degree <- traits$degree[match(spec$PlantGenusSpecies,
                                         traits$GenusSpecies)]

## occurence
occ <- apply(mat, c(3,1), calcOccArray)
spec$occ.date <- apply(spec, 1, findOcc)
traits$occ.date <- spec$occ.date[match(traits$GenusSpecies,
                                       spec$GenusSpecies)]

## plant occurrence
## create sample matrix
site.date <- mat[,,1]
site.date[site.date > 0] <- 0
rownames(site.date) <- lapply(strsplit(rownames(site.date),':'),
                              function(x) x[1])
long.site.date <- comm.mat2sample(site.date)
long.site.date <- long.site.date[!is.na(long.site.date$Samp),]

## create site by date matrices with plant presence
plant.mat <- make.by.species(spec, long.site.date, site.date)
pol.mat <- make.by.species(spec, long.site.date, site.date,
                           type='GenusSpecies')

save(plant.mat, pol.mat, file='../data/species/allSamples.Rdata')

occ.plant <- apply(plant.mat, c(3,1), calcOccArray)

## match to dataset!
spec$occ.plant.date <- apply(spec, 1, findOccPlant)

## update traits file

## mean occurrence accross sites
occ.bees <- apply(occ, 1, mean)

occ.plants <- apply(occ.plant, 1, mean)

occ.all <- c(occ.bees, occ.plants)

traits$occ.date <- occ.all[match(traits$GenusSpecies,
                                 names(occ.all))]

## functional traits
traits <- merge(traits, bee.trait[,c(1:5,27)], all.x=TRUE)

traits <- merge(traits, both.trait[,c(2:3,7)], all.x=TRUE)

traits$bee.syr <- spec$BeeNonbee[match(traits$GenusSpecies,
                                       spec$GenusSpecies)]

traits$bee.syr[is.na(traits$bee.syr)] <- 'plant'

## mean abundance across sites
mean.abund.pol <- apply(apply(pol.mat, c(1,3), mean, na.rm=TRUE), 2,
                        mean)
mean.abund.plant <- apply(apply(plant.mat, c(1,3), mean, na.rm=TRUE),
                          2, mean)
mean.abund <- c(mean.abund.pol, mean.abund.plant)

traits$mean.abun.net <- mean.abund[match(traits$GenusSpecies,
                                         names(mean.abund))]


save(spec, file='../data/networks/allSpecimens.Rdata')
write.csv(traits, file='../data/traits.csv', row.names=FALSE)


site.years <- aggregate(Year~ Site, data=spec,
                        function(x) length(unique(x)))

sites.to.keep <- site.years$Site[site.years$Year >= 1]

## *******************************************************************
## create networks
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
specs <- calcSpec(nets, spec, spec.metric = 'd', 0.3)
specs$closeness[specs$closeness == 0] <- 1*10^-6
specs$closeness.log <- log(specs$closeness)

## within a year, across sites
specs.years <- calcSpec(nets.year, spec, spec.metric = 'd', 0.3)
## within a site across years
specs.site <- calcSpec(nets.site, spec, spec.metric = 'd', 0.3)

specs.save.path <- '../analysis/speciesLevel/saved'
save(specs, specs.years, specs.site,
     file=file.path(specs.save.path, 'specs.Rdata'))



