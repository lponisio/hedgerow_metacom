## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/occupancy')
rm(list=ls())
library(igraph)
library(bipartite)
library(RColorBrewer)
source('plotting/src/makeNetworkFig.R')
source('src/initialize.R')
source('plotting/src/initialize.R')

## only works when keeping track of site level phi and gamma
## for(i in 1:length(phi.gam)){
##     this.site.ave <- phi.gam[[i]]
##     file.name <- c("phi", "gam")[i]

##     pdf.f(plotAllStatuses,
##           file=file.path(save.dir,
##                          sprintf("figures/networks/%s_spatial.pdf",
##                                  file.name)),
##           height=6, width=4)

##     pdf.f(plotbyStatus,
##           file=file.path(save.dir,
##                          sprintf("figures/networks/%s_bystatus.pdf",
##                                  file.name)),
##           height=6, width=8)

##     ## pdf.f(plotbySpecies,
##     ##       file=file.path(save.dir, "figures/networks/byspecies.pdf"),
##     ##       height=6, width=6)
## }


this.site.ave <- NULL
file.name <- "degree"

pdf.f(plotAllStatuses,
      file=file.path(save.dir,
                     sprintf("figures/networks/%s_spatial.pdf",
                             file.name)),
      height=6, width=4)

pdf.f(plotbyStatus,
      file=file.path(save.dir,
                     sprintf("figures/networks/%s_bystatus.pdf",
                             file.name)),
      height=6, width=8)



this.site.ave <- site.between.mean
file.name <- "betweenness"

pdf.f(plotAllStatuses,
      file=file.path(save.dir,
                     sprintf("figures/networks/%s_spatial.pdf",
                             file.name)),
      height=6, width=4)

pdf.f(plotbyStatus,
      file=file.path(save.dir,
                     sprintf("figures/networks/%s_bystatus.pdf",
                             file.name)),
      height=6, width=8)

