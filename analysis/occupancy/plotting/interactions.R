## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom/')
rm(list=ls())
setwd('analysis/occupancy')
library('nimble')
library('RColorBrewer')
source('src/misc.R')
source('src/plotting.R')
source('src/prep.R')
source('src/initialize.R')
source('plotting/src/plotInteractions.R')
source('../../../occupancy/analysis/all/plotting.R')

## *****************************************************************
## interaction plots
## *****************************************************************
natural.decay <- 350

load(file=file.path(save.dir,
                    "runs/mus.Rdata"))
load(file=file.path(save.dir, sprintf("5-0-%s.Rdata", natural.decay)))

means <- mus["mean",]

save(means, file=file.path(save.dir, "means.Rdata"))

probs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.95, 1)
quantiles <- lapply(model.input$data[c("k","B", "HRarea", "natural", "fra")],
                    function(x) quantile(x, probs=probs, na.rm=TRUE))

cols <- brewer.pal(length(probs) + 2, "Blues")[-c(1,2)]

pdf.f(plotInteractions, file=file.path(save.dir,
                                       sprintf("figures/ms/interactions-%s.pdf",
                                               natural.decay)),
      width=9, height=4)


pdf.f(plotHRInteractions, file=file.path(save.dir,
                                         sprintf("figures/ms/HRinteractions-%s.pdf",
                                                 natural.decay)),
      width=9, height=4)



pdf.f(plotInteractionsLinear, file=file.path(save.dir,
                                             sprintf("figures/ms/interactionsLinear-%s.pdf",
                                                     natural.decay)),
      width=9, height=4)


pdf.f(plotAllInteractions, file=file.path(save.dir,
                                          sprintf("figures/ms/all_interactions.pdf",
                                                  natural.decay)),
      width=9, height=8)




