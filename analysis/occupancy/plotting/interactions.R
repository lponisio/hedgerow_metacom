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

load(file=file.path(save.dir,
                    "runs/mus.Rdata"))
load(file=file.path(save.dir, "5-0-all.Rdata"))

means <- mus["mean",]

save(means, file=file.path(save.dir, "means.Rdata"))

probs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.95, 1)
quantiles <- lapply(model.input$data[c("k","B", "HRarea", "natural", "fra")],
                    function(x) quantile(x, probs=probs, na.rm=TRUE))

cols <- brewer.pal(length(probs) + 2, "Blues")[-c(1,2)]

pdf.f(plotInteractions, file=file.path(save.dir,
                                       "figures/ms/interactions.pdf"),
      width=9, height=4)



pdf.f(plotInteractionsLinear, file=file.path(save.dir,
                                       "figures/ms/interactionsLinear.pdf"),
      width=9, height=4)


pdf.f(plotAllInteractions, file=file.path(save.dir,
                                       "figures/ms/all_interactions.pdf"),
      width=9, height=8)




