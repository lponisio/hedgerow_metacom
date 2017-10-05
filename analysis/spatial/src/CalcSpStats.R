library(rgdal)
library(raster)
library(rgeos)

calcSpStats <- function(i,
                        d, ## buffer RADIUS
                        plt, ## plot data
                        nplot, ## number of plots
                        rast, ## raster to calculate stats from
                        ## stats on buffers
                        giv.names=TRUE, ## give the output names
                        plt.name="PLOT" ## name of the column to name
                        ## output with
                        ){
    ## creates buffers of a given size, calculates frag stat metrics
    spStats <-numeric(length(d))
    these.coord <- coordinates(plt[i, ])
    for(j in 1:length(d)){
        p <- spatstat:::disc(d[j], these.coord)
        p <- as(p, 'SpatialPolygons')
        proj4string(p) <- CRS(proj4string(rast))
        ## masking is more time intensive on larger
        ## rasters so crop first
        new.rast <- crop(rast, extent(p))
        if(class(new.rast) == "NULL"){
            spStats[j] <- 0
        } else{
            spStats[j] <- gLength(new.rast)
        }
    }
    return(spStats)
}


applyDecay <- function(dist.tab, decay){
    dist.tab <- dist.tab[[1]]
    areas <- c(dist.tab$area[1], diff(dist.tab$area))
    rads <- c(0, dist.tab$radius[1:(nrow(dist.tab)-1)])
    decay.area <- sum(exp(-(rads-0)^2/(2*decay^2)) * areas)
    return(decay.area)
}
