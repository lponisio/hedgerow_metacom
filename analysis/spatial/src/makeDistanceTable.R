library(parallel)
library(maptools)
library(rgdal)
library(rgeos)
library(spatstat)

makeDistanceTable <- function(dd.lc, dd.h, radii, num.cores,
                              site.col='df0', type.col='VegName',
                              area.function = gArea, drop.self=FALSE) {

    get.areas <- function(ss) {
        cat(sprintf('site: %s\n', ss))
        site.pt <- dd.h[dd.h@data[, site.col]==ss,]

        f.cover.type <- function(type) {
            cat(sprintf('cover type: %s\n', type))
            dd.type <- dd.lc[dd.lc@data[,type.col] == type,]
            if(drop.self){
                dd.type <- dd.type[dd.type@data$site != as.character(ss),]
            }
            ## make list of easily accessible polygons
            ## make.poly <- function(poly){
            ##     SpatialPolygons(list(Polygons(list(poly), poly@ringDir)),
            ##                     proj4string=CRS(proj4string(dd.type)))
            ## }
            ## polys.list <- sapply(dd.type@polygons[[1]]@Polygons, make.poly)

            ## drop tiny polygons
            ## polys.list <- polys.list[sapply(polys.list, area) > 1]

            ## add zero width buffers around polygons (this fixes problems
            ## below)
            ## polys.list <- sapply(polys.list,
            ##                      function(x) gBuffer(x, width=0, byid=TRUE))

            get.intersection.area <- function(radius) {
                ## create buffer around site
                buff <- gBuffer(site.pt, width=radius)
                new.poly <- crop(dd.type, buff)
                if(is.null(new.poly)){
                    area.nat <- 0
                } else{
                    area.nat <- area.function(new.poly)
                }

                ## sum.areas <- ifelse(length(intersections)==0,
                ##                     0, sum(sapply(intersections, area)))
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

    dist.tables
}
