library(parallel)
library(maptools)
library(rgdal)
library(rgeos)
library(spatstat)

makeDistanceTable <- function(dd.lc, dd.h, radii, num.cores,
                              site.col='df0', type.col='VegName',
                              area.function = gArea, drop.self=FALSE) {
    ## dd.lc = data to calculate area from
    ## dd.h = focal site points
    ## radii = vector of radii to create buffers to calculate area
    ## within
    ## site.col = site column withitn @data with site names
    ## type.col = type column within dd.lc to calculate area of
    ## area.function = gArea for natural area, and gLength for
    ## hedgerows
    ## drop.self, drop sites from the data area is being calculated on

    get.areas <- function(ss) {
        cat(sprintf('site: %s\n', ss))
        site.pt <- dd.h[dd.h@data[, site.col]==ss,]

        f.cover.type <- function(type) {
            cat(sprintf('cover type: %s\n', type))
            dd.type <- dd.lc[dd.lc@data[,type.col] == type,]
            ## drop focal site from data used to calculate area
            if(drop.self){
                dd.type <- dd.type[dd.type@data$site != as.character(ss),]
            }

            get.intersection.area <- function(radius) {
                ## create buffer around site
                buff <- gBuffer(site.pt, width=radius)
                new.poly <- crop(dd.type, buff)
                if(is.null(new.poly)){
                    ## if no area within buffer
                    area.nat <- 0
                } else{
                    area.nat <- area.function(new.poly)
                }
                data.frame(area=area.nat,
                           site=ss,
                           radius=radius)
            }
            do.call(rbind, lapply(radii, get.intersection.area))
        }
        res <- lapply(unique(dd.lc@data[, type.col]), f.cover.type)
        names(res) <- unique(dd.lc@data[, type.col])
        res
    }

    ## create list of areas for each site
    site.ids <- dd.h@data[, site.col]
    if(num.cores==1)
        dist.tables <- lapply(site.ids, get.areas)
    if(num.cores > 1)
        dist.tables <- mclapply(site.ids, get.areas,
                                mc.preschedule=FALSE,
                                mc.cores=num.cores)
    names(dist.tables) <- site.ids

    return(dist.tables)
}
