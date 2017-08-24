

## calculate spatial type area buffers, create natural data

calcWeightedCover <- function(yr, site, type, decay, dist.tables) {
    dist.tab <- dist.tables[[as.character(yr)]][[site]][[type]]
    areas <- c(dist.tab$area[1], diff(dist.tab$area))
    rads <- c(0,dist.tab$radius[1:(nrow(dist.tab)-1)])
    sum(exp(-(rads-0)^2/(2*decay^2)) * areas)
}

## decay 350, 1000, 2500
getMetric <- function(vv, dd, gis.yr, decay, dist.tables)
    mapply(function(a, b)
        calcWeightedCover(yr=a, site=b, type=vv,
                            decay=decay,
                            dist.tables=dist.tables),
        a=gis.yr, b=dd$Site)

getNatural <- function(spec, dist.tables, decays){
    cols <-  paste0(c("riparianForest", "riparianScrub", "natural"),
                    rep(decays, each=3))

    natural <- unique(data.frame(Site=spec$Site,
                                 Year=spec$Year))
    natural <- cbind(natural, matrix(NA, ncol=length(cols),
                                     nrow=nrow(natural)))
    colnames(natural) <- c(colnames(natural)[1:2], cols)
    gis.yr <- natural$Year
    gis.yr[gis.yr < 2010] <- 2006
    gis.yr[gis.yr >= 2010] <- 2012

    for(decay in decays){
        natural[, paste0("riparianForest", decay)] <-
            getMetric(vv='Riparian Forest',
                       dd=natural,
                       gis.yr=gis.yr,
                       decay=decay,
                       dist.tables =dist.tables)
        natural[, paste0("riparianScrub", decay)] <-
            getMetric(vv='Riparian Scrub',
                       dd=natural,
                       gis.yr=gis.yr,
                       decay=decay,
                       dist.tables =dist.tables)
        natural[, paste0("natural", decay)] <- natural[,
                                      paste0("riparianForest", decay)]+
            natural[, paste0("riparianScrub", decay)]
    }
    return(natural)
}
