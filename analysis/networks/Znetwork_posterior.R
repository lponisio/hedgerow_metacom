## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')
rm(list=ls())
source('src/initialize.R')
source('src/zNetPrep.R')
library(parallel)
options(ncores=10)

## load the GIANT model file that tracks all of the latent states
load.dir <- "../../../hedgerow_metacom_saved/occupancy/runs"
load(file.path(load.dir, 'all_2500_350_FALSE.Rdata'))

zs <- grepl("Z", colnames(ms.ms.nimble[[1]]))
samples.z <- lapply(ms.ms.nimble, function(x) x[,zs])

## take a random sample of samples from a random chain
these.samples <- samples.z[[sample(1:3, 1)]]
these.samples <- these.samples[sample(1:nrow(these.samples),
                                      size=1000, replace=FALSE),]

year.nets <- mclapply(1:nrow(these.samples), calcYearNet,
                      these.samples, model.input)

## metacommunity patch network
year.site <- mclapply(1:dim(model.input$data$Z)[2], getPosteriorNet,
                               lapply(year.nets, function(x) x$site),
                    model.input, "site")
## pollinator species network
year.sp <- mclapply(1:dim(model.input$data$Z)[2], getPosteriorNet,
                               lapply(year.nets, function(x) x$sp),
                    model.input, "sp")
## species temporal network
site.sp <- lapply(1:dim(model.input$data$Z)[1], getPosteriorNet,
                               lapply(year.nets, function(x) x$year),
                    model.input, "year")

## add trait data
year.site <- do.call(rbind, year.site)
year.sp <- do.call(rbind, year.sp)
site.sp <- do.call(rbind, site.sp)
year.sp <-  merge(year.sp,
                       traits[,c("GenusSpecies", "r.degree",
                                 "MeanITD")])
site.sp <-  merge(site.sp,
                       traits[,c("GenusSpecies", "r.degree",
                                 "MeanITD")])
year.site <- merge(year.site, by.site)

## **********************************************************
## pollinators
## **********************************************************

## anything outputted by specieslevel
ys <- c("degree", "betweenness", "closeness")
mod.weights <- c("sd.degree", "sd.betweenness", "sd.closeness")

xvar.species.form <- paste0("scale(", xvar.species, ")")

formulas.species <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.species.form , collapse="+"),
                           "(1|GenusSpecies)",
                           sep="+")))
})


## within a year, across sites
mod.years.pol <- list()
mod.years.pol[[1]] <- lmer(formulas.species[[1]], data=year.sp,
                               weights=sqrt((1/year.sp$sd.degree)))

mod.years.pol[[2]] <- lmer(formulas.species[[2]], data=year.sp,
                           weights=sqrt((1/year.sp$sd.betweenness)))

mod.years.pol[[3]] <- lmer(formulas.species[[3]], data=year.sp,
                           weights=sqrt((1/year.sp$sd.closeness)))

## within a site across years
mod.sites.pol <- list()
mod.sites.pol[[1]] <- lmer(formulas.species[[1]], data=site.sp,
                           weights=sqrt((1/site.sp$sd.degree)))

mod.sites.pol[[2]] <- lmer(formulas.species[[2]], data=site.sp,
                           weights=sqrt((1/site.sp$sd.betweenness)))

mod.sites.pol[[3]] <- lmer(formulas.species[[3]], data=site.sp,
                           weights=sqrt((1/site.sp$sd.closeness)))


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
xvar.site.form <- paste0("scale(", xvar.site, ")")
formulas.site <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.site.form, collapse="+"),
                           "(1|Site)",
                           sep="+")))
})
mod.years.site <- list()
mod.years.site[[1]] <- lmer(formulas.site[[1]], data=year.site,
                           weights=sqrt((1/year.site$sd.degree)))

mod.years.site[[2]] <- lmer(formulas.site[[2]], data=year.site,
                            weights=sqrt((1/year.site$sd.betweenness)))

mod.years.site[[3]] <- lmer(formulas.site[[3]], data=year.site,
                           weights=sqrt((1/year.site$sd.closeness)))


## name them the same as the pollinators
names(mod.years.site)  <- ys

print("metacommunity patch network")
print(lapply(mod.years.site, summary))

save(mod.years.pol, mod.sites.pol, mod.years.site,
     year.sp, site.sp, year.site,
     file=sprintf('saved/mods/zmets_drop_li_ht%s.Rdata', drop.li.ht))
