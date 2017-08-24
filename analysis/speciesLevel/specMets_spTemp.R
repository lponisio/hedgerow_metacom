rm(list=ls())
## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/speciesLevel')
xvar <- c("r.degree", "BodyLength")
source('src/initialize.R')

## **********************************************************
## species level metrics
## **********************************************************


## anything outputted by specieslevel
ys <- c("proportional.generality.x",  "degree.x", "k")

formulas <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar, collapse="+"),
                           "(1|Site)",
                           "(1|GenusSpecies)",
                           sep="+")))
})

mod.sites <- lapply(formulas, function(x){
    lmer(x, data=specs.site.pol)
})


mod.years <- lapply(formulas, function(x){
    lmer(x, data=specs.years.pol)
})


names(mod.years) <- names(mod.sites) <- ys

lapply(mod.years, summary)
lapply(mod.sites, summary)

save(mod.years, mod.sites,
     file='saved/mods/specmetrics.Rdata')
