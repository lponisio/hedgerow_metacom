## ************************************************************
rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/occupancy')
library('abind')
library('nimble')
library('R2jags')
source('src/misc.R')
source('src/plotting.R')
source('src/prep.R')
source('src/initialize.R')
source('../../../occupancy/analysis/all/plotting.R')

## *****************************************************************
## comparisons
## *****************************************************************
## load(file=file.path(save.dir, 'runs/crosslevel.Rdata'))
load(file=file.path(save.dir, 'runs/nimble_bees_noRain.Rdata'))


ms.ms.occ.all <- combine_MCMC_comparison_results(ms.ms.nimble[[1]],
                                                 ## ms.ms.crosslevel[[1]],
                                                 name = "ms.ms")

## doesn't really work with a lot of params
## make_MCMC_comparison_pages(ms.ms.occ.all,
##                            dir=file.path(save.dir, "figures/comparisons"))

## takes forever with a lot of params
## checkChains(ms.ms.occ.all$ms.ms$samples,
##             f.path = file.path(save.dir, "figures/chains/%s.pdf"))

## *****************************************************************
library(ggmcmc)

## mus
samples.mus <- ggs(as.mcmc(t(ms.ms.occ.all$ms.ms$samples["nimble",
  grep("^mu", dimnames(ms.ms.occ.all$ms.ms$samples["nimble",,])[[1]]),])))

ggmcmc(samples.mus, file=file.path(save.dir,
                                   "figures/chains/nimble-mus.pdf"))


cat.func <- function() {plot(ggs_caterpillar(samples.mus))}

pdf.f(cat.func,
      file=file.path(save.dir,
                     "figures/chains/nimble-mus-caterpillar.pdf"),
      height=6, width=4)

## *****************************************************************
## paramter groups
## *****************************************************************

params <- dimnames(ms.ms.occ.all$ms.ms$summary)[[3]]
groups <- list()
i <- 1
while(length(params) > 0){
    id <- agrep(params[1], params, max.distance = 0.3)
    groups[[i]] <- params[id]
    params <- params[-id]
    i <- i + 1
}

pdf.f(plotComparisons,
      file=file.path(save.dir,
                     "figures/comparisons/allparams.pdf"),
      height=10, width=8.5 )


## *****************************************************************
## posterior plots for ms
## *****************************************************************

nimble.summary <- ms.ms.occ.all$ms.ms$summary["nimble",,]

mus <- nimble.summary[,grep("^mu", colnames(nimble.summary))]

wanted.order <- c("hr.area", "nat.area", "fra", "hr.area.fra",
                  "nat.area.fra", "traits1", "traits1.fra", "traits2")

xlabs <- c("Hedgerow proximity", "Semi-natural \n habitat proximity",
           "Floral diversity", "Hedgerow* \n floral diversity",
           "Semi-natural* \n floral diversity", "Diet breadth",
           "Diet breadth* \n floral diversity",
           "Body size")

f <- function() {plotPosterior(mus, wanted.order, xlabs)}

pdf.f(f,
      file=file.path(save.dir,
                     "figures/ms/mus.pdf"),
      height=8, width=6)



## *****************************************************************
## persistence vs. hr and nat habitat effects
## *****************************************************************
load(file.path(save.dir,
               'runs/nimble_bees_noRain_allparam_short.Rdata'))

load(file.path(save.dir,
               '5-0-all.Rdata'))

nimble.sum <- ms.ms.nimble$model1$summary["nimble",,]

params.to.get <- c("phi.sp.mean", "gam.sp.mean", "phi.nat.area",
                   "gam.nat.area", "phi.hr.area", "gam.hr.area")

params <- lapply(params.to.get, function(x){
    nimble.sum[, grep(x,
                      colnames(nimble.sum))]
})
names(params) <- params.to.get

params <- lapply(params, function(x){
    x[, !grepl("^mu|^sigma", colnames(x))]
})

## no relationships between persistence/colonization and hedgerow/nat
## habtiat proximity effect
layout(matrix(1:4, byrow=TRUE, nrow=2))
plot(params$phi.sp.mean["mean",], params$phi.hr.area["mean",])
plot(params$gam.sp.mean["mean",], params$gam.hr.area["mean",])

plot(params$phi.sp.mean["mean",], params$phi.nat.area["mean",])
plot(params$gam.sp.mean["mean",], params$gam.nat.area["mean",])

