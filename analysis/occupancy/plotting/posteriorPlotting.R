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
                                          natural.decay,
                            include.int)))
load(file=file.path(save.dir, sprintf("5-0-%s.Rdata", natural.decay)))


if(is.list(ms.ms.nimble$samples)){
    all.samples <- do.call(rbind, ms.ms.nimble$samples)
} else{
    all.samples <- ms.ms.nimble$samples
}

nimble.summary <- apply(all.samples, 2, function(x){
    means <- mean(x)
    CI95_upp <- 1.96*sd(x) + means
    CI95_low <-  means - 1.96*sd(x)
    return(c(mean=means, CI95_upp=CI95_upp, CI95_low=CI95_low ))
})

mus <- nimble.summary[,grep("^mu", colnames(nimble.summary))]

wanted.order <- c("hr.area",
                  "nat.area",
                  "fra",
                  "k",
                  "B",
                  "hr.area.fra",
                  "nat.area.fra",
                  "hr.area.k",
                  "nat.area.k",
                  "hr.area.B",
                  "nat.area.B")

xlabs <- c("Hedgerow \n area/proximity",
           "Non-crop habitat \n area/proximity",
           "Floral diversity",
           "Floral diet breadth",
           "Body size",
           "Hedgerow \n area/proximity*\n floral diversity",
           "Non-crop \n area/proximity*\n floral diversity",
           "Hedgerow \n area/proximity*\n floral diet breadth",
           "Non-crop \n area/proximity*\n floral diet breadth",
           "Hedgerow \n area/proximity*\n body size",
           "Non-crop \n area/proximity*\n body size")


if(include.int == "no_noncrop"){
    wanted.order <- wanted.order[!grepl("nat", wanted.order)]
    xlabs <- xlabs[!grepl("Non-crop", xlabs)]
}

## variables to plot
by.site <- by.site[by.site$Site %in% rownames(model.input$data$fra),]
controls <- by.site$Site[by.site$SiteStatus == "control"]
hedgerows <- by.site$Site[by.site$SiteStatus == "mature" | by.site$SiteStatus == "maturing"]


f <- function() {plotPosterior(mus, wanted.order, xlabs)}

pdf.f(f,
      file=file.path(save.dir,
                       sprintf('figures/posterior/mus_bees_%s_%s.pdf',
                                          natural.decay, include.int)),
      height=7, width=9)

save(mus, file=file.path(save.dir,
                      sprintf('runs/mus_bees_%s_%s.Rdata',
                                          natural.decay, include.int)))


all.samples <- all.samples[, colnames(all.samples) %in%
                                   ms.ms.model$getNodeNames(includeData=FALSE,
                                                            stochOnly=TRUE)]



## runMCMCcheckChains(ms.ms.nimble$samples, f.path=file.path(save.dir,
##                                                   "figures/chains"),
##                                                   natural.decay, include.int)

## ****************************************************************
##
## posterior probability table
##
## *****************************************************************

load(file=file.path(save.dir,
                    sprintf('runs/post_probs_nimble_bees_%s_%s.Rdata',
                            natural.decay, include.int)))

phis <- paste("mu.phi", wanted.order,
              sep=".")
gams <- paste("mu.gam", wanted.order,
              sep=".")

probs.4.table.phis <- round(posterior.probs[rownames(posterior.probs) %in%
                                      phis,],3)[phis,]
probs.4.table.gams <- round(posterior.probs[rownames(posterior.probs) %in%
                                      gams,], 3)[gams,]

rownames(probs.4.table.phis) <- rownames(probs.4.table.gams) <- xlabs

write.table(cbind(probs.4.table.phis[,-2], probs.4.table.gams[,-2]),
            sep=" & ",
            file=file.path(save.dir, "table/posteriorProbs.txt"))
