## ************************************************************
## setwd('~/Dropbox/hedgerow_metacom/')
rm(list=ls())
setwd('analysis/occupancy')
library('nimble')
source('src/misc.R')
source('plotting/src/plotting.R')
source('plotting/src/checkChains.R')
source('src/prep.R')
source('src/initialize.R')
include.int <- "allInt"
## include.int <- "no_noncrop"
## 350, 1000, 2500
natural.decay <- "350"

load(file=file.path(save.dir,
                    sprintf('runs/nimble_bees_%s_%s.Rdata',
                                          natural.decay, include.int)))

all.samples <- do.call(rbind, ms.ms.nimble$samples)

nimble.summary <- apply(all.samples, 2, function(x){
    means <- mean(x)
    CI95_upp <- 1.96*sd(x) + means
    CI95_low <-  means - 1.96*sd(x)
    return(c(mean=means, CI95_upp=CI95_upp, CI95_low=CI95_low ))
})

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


if(include.int == "no_noncrop"){
    wanted.order <- wanted.order[!grepl("nat", wanted.order)]
    xlabs <- xlabs[!grepl("Non-crop", xlabs)]
}


f <- function() {plotPosterior(mus, wanted.order, xlabs)}

pdf.f(f,
      file=file.path(save.dir,
                       sprintf('figures/posterior/mus_bees_%s_%s.pdf',
                                          natural.decay, include.int)),
      height=7, width=9)

save(mus, file=file.path(save.dir,
                      sprintf('runs/mus_bees_%s_%s.Rdata',
                                          natural.decay, include.int)))


runMCMCcheckChains(ms.ms.nimble$samples, f.path=file.path(save.dir,
                                                  "figures/chains"),
                                                  natural.decay, include.int)
