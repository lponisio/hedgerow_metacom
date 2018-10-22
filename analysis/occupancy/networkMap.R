## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/occupancy')
rm(list=ls())
source('src_plotting/makeNetworkFig.R')
source('src_plotting/initialize.R')

## this doesn't work sometimes. No idea why.
sys <- gmap(x=bbox.all,
            scale=2,
            type="satellite")

all.sites.pt <- spTransform(all.sites.pt,
                            sys@crs)
lat.long <- coordinates(all.sites.pt)
rownames(lat.long) <- all.sites.pt@data$df0
nets.year.sp <- lapply(nets.year, t)
landcover.nat <- spTransform(landcover.nat,
                             sys@crs)

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



## this.site.ave <- site.between.mean
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

