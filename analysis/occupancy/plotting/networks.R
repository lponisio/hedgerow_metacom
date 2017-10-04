rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/occupancy')
library(igraph)
library(bipartite)
library(RColorBrewer)
source('src/initialize.R')
source('plotting/src/initialize.R')

for(i in 1:length(phi.gam)){
    this.site.ave <- phi.gam[[i]]
    file.name <- c("phi", "gam")[i]

    pdf.f(plotAllStatuses,
          file=file.path(save.dir,
                         sprintf("figures/networks/%s_spatial.pdf", file.name)),
          height=6, width=4)

    pdf.f(plotbyStatus,
          file=file.path(save.dir,
                         sprintf("figures/networks/%s_bystatus.pdf", file.name)),
          height=6, width=8)

    ## pdf.f(plotbySpecies,
    ##       file=file.path(save.dir, "figures/networks/byspecies.pdf"),
    ##       height=6, width=6)
}
