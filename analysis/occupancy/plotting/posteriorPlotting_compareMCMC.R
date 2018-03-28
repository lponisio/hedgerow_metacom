## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom/')
rm(list=ls())
setwd('analysis/occupancy')
library('nimble')
source('src/misc.R')
source('src/plotting.R')
source('src/prep.R')
source('src/initialize.R')
source('../../../occupancy/analysis/all/plotting.R')
include.int <- "allInt"
## 350, 1000, 2500
natural.decay <- "1000"
data.subset <- "all"

load(file=file.path(save.dir,
                    sprintf('runs/%s_nimble_bees_%s_%s.Rdata',
                                          data.subset,
                                          natural.decay,
                            include.int)))
## *****************************************************************
## comparisons
## *****************************************************************
## if using compareMCMC

ms.ms.occ.all <- combine_MCMC_comparison_results(ms.ms.nimble,
                                                 name = "ms.ms")

## doesn't really work with a lot of params
## make_MCMC_comparison_pages(ms.ms.occ.all,
##                            dir=file.path(save.dir, "figures/comparisons"))

## takes forever with a lot of params
checkChains(ms.ms.occ.all$ms.ms$samples, only.one="nimble",
            f.path = file.path(save.dir,
                                  sprintf('figures/%s_nimble_bees_%s_%s.Rdata',
                                          data.subset,
                                          natural.decay,
                            include.int)))

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

wanted.order <- c("hr.area", "nat.area", "fra", "k", "B", "hr.area.fra",
                  "nat.area.fra",
                  "hr.area.k", "nat.area.k",
                  "hr.area.B", "nat.area.B")

xlabs <- c("Hedgerow proximity", "Non-crop \n habitat proximity",
           "Floral diversity", "Floral diet breadth",
           "Body size",
           "Hedgerow proximity* \n floral diversity",
           "Non-crop proximity* \n floral diversity",
           "Hedgrow proximity* \n floral diet breadth",
           "Non-crop proximity* \n floral diet breadth",
           "Hedgrow proximity* \n body size",
           "Non-crop proximity* \n body size")

f <- function() {plotPosterior(mus, wanted.order, xlabs)}

pdf.f(f,
      file=file.path(save.dir,
                     "figures/ms/mus.pdf"),
      height=7, width=9)

## calculate some thigns for the ms
mus["mean","mu.phi.nat.area.fra"] + mus["mean","mu.phi.nat.area"]


save(mus, file=file.path(save.dir,
                     "runs/mus.Rdata"))

## *****************************************************************
## persistence vs. hr and nat habitat effects
## *****************************************************************
## load(file=file.path(save.dir,
##                     'runs/nimble_bees_noRain_short.Rdata'))

## load(file.path(save.dir,
##                '5-0-all.Rdata'))

## nimble.sum <- ms.ms.nimble$model1$summary["nimble",,]

## params.to.get <- c("phi.sp.mean", "gam.sp.mean",
##                    'phi.nat.area',
##                    'phi.hr.area',
##                    'phi.hr.area.fra',
##                    'phi.nat.area.fra',
##                    'phi.hr.area.k')

## params <- lapply(params.to.get, function(x){
##     nimble.sum[, grep(x,
##                       colnames(nimble.sum))]
## })
## names(params) <- params.to.get

## params <- lapply(params, function(x){
##     x[, !grepl("^mu|^sigma", colnames(x))]
## })

