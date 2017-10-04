rm(list=ls())
## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/speciesLevel')

xvar.species <- c("r.degree", "BodyLength")

xvar.site <- c("Div", "natArea", "hrArea")

natural.decay <- "350" ## match with occupancy model

source('src/initialize.R')

## **********************************************************
## pollinators
## **********************************************************

## anything outputted by specieslevel
ys <- c("proportional.generality.x",  "degree.x", "k", "betweenness", "closeness")


formulas.species <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.species, collapse="+"),
                           "(1|Site)",
                           "(1|GenusSpecies)",
                           sep="+")))
})

## within a site across years
mod.sites.pol <- lapply(formulas.species, function(x){
    lmer(x, data=specs.site.pol)
})

## within a year, across sites
mod.years.pol <- lapply(formulas.species, function(x){
    lmer(x, data=specs.years.pol)
})


names(mod.years.pol) <- names(mod.sites.pol) <- ys

## floral degree positivly related to all measures of network
## importance in space/time
lapply(mod.years.pol, summary)
lapply(mod.sites.pol, summary)

## **********************************************************
## sites
## **********************************************************

formulas.site <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.site, collapse="+"),
                           "(1|Site)",
                           "(1|GenusSpecies)",
                           sep="+")))
})

## within a site across years. This is not really that
## interesting. Why is it imporant to know what year is most central?
mod.sites.site <- lapply(formulas.site, function(x){
    lmer(x, data=specs.site.site)
})

## within a year across sites
mod.years.site <- lapply(formulas.site, function(x){
    lmer(x, data=specs.years.site)
})

## name them the same as the pollinators
names(mod.years.site) <- names(mod.sites.site) <- ys

## natural area and floral div are sig for the degree-like measures of
## importance, but not the centrality

lapply(mod.years.site, summary)
lapply(mod.sites.site, summary)


save(mod.years.pol, mod.sites.pol, mod.years.site, mod.sites.site,
     file='saved/mods/specmetrics.Rdata')
