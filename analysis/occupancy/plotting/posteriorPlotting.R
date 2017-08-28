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

## *****************************************************************
## comparisons
## *****************************************************************
source('../../../occupancy/analysis/all/plotting.R')

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


## phis
samples.phis <- ggs(as.mcmc(t(ms.ms.occ.all$ms.ms$samples["nimble",
 grep("^phi", dimnames(ms.ms.occ.all$ms.ms$samples["nimble",,])[[1]]),])))

ggmcmc(samples.phis,
       file=file.path(save.dir,
                      "figures/chains/nimble-phis.pdf"))

pdf.f(ggs_caterpillar(samples.phis),
      file=file.path(save.dir,
                     "figures/chains/nimble-phis-caterpillar.pdf"))



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
                  "nat.area.fra", "traits1", "traits2") # "traits1.fra"

xlabs <- c("Hedgerow proximity", "Semi-natural \n habitat proximity",
           "Floral diversity", "Hedgerow* \n floral diversity",
           "Semi-natural* \n floral diversity", "Diet breadth",
           "Body size")

f <- function() {plotPosterior(mus, wanted.order, xlabs)}

pdf.f(f,
      file=file.path(save.dir,
                     "figures/ms/mus.pdf"),
      height=8, width=6)

