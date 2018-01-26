## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom/')
rm(list=ls())
setwd('analysis/occupancy')
library('nimble')
library('viridis')
source('src/misc.R')
source('plotting/src/plotting.R')
source('src/prep.R')
source('src/initialize.R')
source('plotting/src/plotInteractions.R')
source('../../../occupancy/analysis/all/plotting.R')

## *****************************************************************
## interaction plots
## *****************************************************************

include.int <- "allInt"
## include.int <- "no_noncrop"
## 350, 1000, 2500
natural.decay <- "350"

load(file=file.path(save.dir,
                      sprintf('runs/mus_bees_%s_%s.Rdata',
                                          natural.decay, include.int)))


load(file=file.path(save.dir, sprintf("5-0-%s.Rdata", natural.decay)))

means <- mus["mean",]
inv.logit(means)

save(means, file=file.path(save.dir, "means.Rdata"))
probs <- seq(from=0, to=1, by=0.05)
quantiles <- lapply(model.input$data[c("k","B", "HRarea", "natural", "fra")],
                    function(x) quantile(x, probs=probs, na.rm=TRUE))

pdf.f(plotHRPersistence, file=file.path(save.dir,
                                       sprintf("figures/interactions/HRpersistence-%s.pdf",
                                               natural.decay)),
      width=9, height=4)

pdf.f(plotHRPerstTurnOccCol, file=file.path(save.dir,
                                         sprintf("figures/interactions/HRinteractions-%s.pdf",
                                                 natural.decay)),
      width=9, height=8)
