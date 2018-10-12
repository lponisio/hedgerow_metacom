## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')
source('src/initialize.R')

## **********************************************************
## pollinators
## **********************************************************

## anything outputted by specieslevel
ys <- c("k", "weighted.betweenness")
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
print("metacommunity pollinator spatial network")
print(lapply(mod.years.pol, summary))
print("metacommunity pollinator temporal network")
print(lapply(mod.sites.pol, summary))

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

## within a year across sites
mod.years.site <- lapply(formulas.site, function(x){
    lmer(x, data=specs.years.site)
})

## name them the same as the pollinators
names(mod.years.site)  <- ys

print("metacommunity patch network")
print(lapply(mod.years.site, summary))

save(mod.years.pol, mod.sites.pol, mod.years.site,
     file=sprintf('saved/mods/specmetrics_drop_li_ht%s.Rdata', drop.li.ht))
