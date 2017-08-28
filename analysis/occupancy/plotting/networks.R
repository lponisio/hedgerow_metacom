rm(list=ls())
setwd('~/Dropbox/hedgerow_metacom/analysis/occupancy')
library(igraph)
library(bipartite)
library(RColorBrewer)
source('src/initialize.R')
source('src/misc.R')
source('plotting/src/makeNetworkFig.R')



plotAllStatuses <- function(){
    makeNetworkFig(spec, sys, lat.long=lat.long,
                   rows="Site",
                   cols="GenusSpecies",
                   nets=nets.year,
                   natural.cover=landcover.nat,
                   add=TRUE,
                   rescale = FALSE)
}



plotbySpecies <- function(){
    makeNetworkFig(spec,
                   rows="GenusSpecies",
                   cols="Site",
                   nets=nets.year.sp,
                   natural.cover=landcover.nat,
                   layout=layout_nicely,
                   edge.color=add.alpha("grey", 0.1),
                   rescale=TRUE)
}


plotbyStatus <- function(){
    layout(matrix(1:3, ncol=3))
    cols.e <- add.alpha(brewer.pal(11, 'RdGy')[c(6,11,1,6)],
                        alpha=0.5)
    cols.v <- brewer.pal(11, 'RdGy')[c(1,11)]
    print("mature")
    makeNetworkFig(spec[spec$SiteStatus == "mature",],
                   sys, lat.long=lat.long,
                   rows="Site",
                   cols="GenusSpecies",
                   nets=nets.year,
                   cols.vertex=cols.v[2],
                   cols.edges=cols.e[2],
                   vertex.names="mature",
                   natural.cover=landcover.nat,
                   add=TRUE,
                   rescale = FALSE)
    print("control")
    makeNetworkFig(spec[spec$SiteStatus == "control",],
                   sys, lat.long=lat.long,
                   rows="Site",
                   cols="GenusSpecies",
                   nets=nets.year,
                   cols.vertex=cols.v[1],
                   cols.edges=cols.e[3],
                   vertex.names="control",
                   natural.cover=landcover.nat,
                   add=TRUE,
                   rescale = FALSE)

    cols.e[c(2:3)] <- add.alpha("white", 0)
    print("mixed")
    makeNetworkFig(spec,
                   sys, lat.long=lat.long,
                   rows="Site",
                   cols="GenusSpecies",
                   nets=nets.year,
                   natural.cover=landcover.nat,
                   cols.vertex=cols.v,
                   cols.edges=cols.e,
                   add=TRUE,
                   rescale = FALSE)
}


pdf.f(plotAllStatuses,
      file=file.path(save.dir, "figures/networks/spatial.pdf"),
      height=6, width=4)

pdf.f(plotbyStatus,
      file=file.path(save.dir, "figures/networks/bystatus.pdf"),
      height=6, width=8)

pdf.f(plotbySpecies,
      file=file.path(save.dir, "figures/networks/byspecies.pdf"),
      height=6, width=6)
