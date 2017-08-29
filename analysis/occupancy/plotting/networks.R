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
    plotAllStatuses <- function(){
        makeNetworkFig(spec, sys, lat.long=lat.long,
                       rows="Site",
                       cols="GenusSpecies",
                       nets=nets.year,
                       natural.cover=landcover.nat,
                       add=TRUE,
                       rescale = FALSE,
                       site.ave=this.site.ave)
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
        cols.e <- add.alpha(brewer.pal(11, 'RdGy')[c(6,10,1,6)],
                            alpha=0.5)
        cols.v <- brewer.pal(11, 'RdGy')[c(1,10)]
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
                       rescale = FALSE,
                       site.ave=this.site.ave)
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
                       rescale = FALSE,
                       site.ave=this.site.ave)

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
                       rescale = FALSE,
                       site.ave=this.site.ave)
    }


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
