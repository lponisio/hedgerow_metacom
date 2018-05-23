## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom/')
rm(list=ls())
setwd('analysis/occupancy')
library('nimble')
library('viridis')
source('src/misc.R')
source('plotting/src/plotting.R')
source('src/prep.R')
args <- commandArgs(trailingOnly=TRUE)
args <- c("allInt","350","filtering","all",2e2)
source('src/initialize.R')
source('plotting/src/plotInteractions.R')
source('../../../occupancy/analysis/all/plotting.R')

## *****************************************************************
## interaction plots
## *****************************************************************

load(file=file.path(save.dir,
                        sprintf('runs/%s_nimble_bees_%s_%s.Rdata',
                                data.subset, natural.decay, include.int)))

load(file=file.path(save.dir, sprintf("%s-5-0-%s.Rdata", data.subset,
                                      natural.decay)))

load(file=file.path(save.dir,
                             sprintf('runs/%s_mus_bees_%s_%s.Rdata',
                                     data.subset, natural.decay, include.int)))


if(is.list(ms.ms.nimble$samples)){
    all.samples <- do.call(rbind, ms.ms.nimble$samples)
} else{
    all.samples <- ms.ms.nimble$samples
}

means <- apply(all.samples, 2, mean)

pdf.f(plotInteractionsFloralDiv, file=file.path(save.dir,
                                         sprintf("figures/interactions/HRinteractions-fra-%s.pdf",
                                                 natural.decay)),
      width=5, height=9)



pdf.f(plotInteractionsK, file=file.path(save.dir,
                                         sprintf("figures/interactions/HRinteractions-k-%s.pdf",
                                                 natural.decay)),
      width=4, height=11)



pdf.f(plotInteractionsB, file=file.path(save.dir,
                                         sprintf("figures/interactions/HRinteractions-B-%s.pdf",
                                                 natural.decay)),
      width=4, height=11)




se <- apply(all.samples, 2, sd)
