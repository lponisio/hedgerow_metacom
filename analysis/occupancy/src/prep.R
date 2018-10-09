
prepOccModelInput <- function(nzero, ## if augmenting data
                              threshold, ## min times a pollinator is seen
                              spec, ## specimen data
                              sr.sched, ## sampling schedule
                              traits, ## dataframe with traits
                              col.name.trait1, ## colname of one trait
                              col.name.trait2, ## colname of second trait
                              HRarea, ## vector of hedgerow area prox
                              natural.mat=NULL, ## site by year
                              ## natural area proximity
                              natural.decay,## site by decay matrix
                              ## natural
                              HR.decay,## site by decay matrix natura
                              veg, ## floral availability matrix
                              col.name.div.type="Div", ## colname of
                              ## the metric of floral availability
                              use.HR.decay=TRUE,
                              data.subset="all",
                              save.dir) {

    ## this function preps specimen data for the multi season multi
    ## species occupancy model, and returns a lsit of the inits, data,
    ## constants needed to run in jags or nimble

    ## create complete site x date x species matrix based on sample
    ## schedule
    null.mat <- tapply(rep(0, nrow(sr.sched)),
                       list(sites=paste(sr.sched$Site),
                            dates=sr.sched$Date), sum)

    ## extract the relevant specimen, site, date data
    spec.data <- data.frame(pollinator=spec$GenusSpecies,
                            site=spec$Site,
                            date=spec$Date)
    ## natural cover
    natural.mat <- natural.mat[, natural.decay]
    natural.mat <- natural.mat[!is.na(natural.mat)]
    natural.sites <- names(natural.mat)

    ## area of hedgerow
    if(use.HR.decay){
        hr.mat <- HRarea[, HR.decay]
        d.area <- hr.mat[!is.na(hr.mat)]
        buffer <- "all"
    } else{
        d.area <- HRarea
        buffer <- "all"
    }
    ## log and standardize HR area
    d.area <- log(d.area[names(d.area) %in%
                         natural.sites] + 1)
    d.area <- standardize(d.area)

    flower.mat <- samp2site.ypr(site=veg$Site,
                                yr=veg$Year,
                                abund=veg[, col.name.div.type])

    flower.mat <- t(apply(flower.mat, 1, function(x){
        nas <- which(is.na(x))
        if(length(nas) != 0){
            x[is.na(x)] <- mean(x, na.rm=TRUE)
        }
        return(x)
    }))

    ## drop sites without natural data or hedgerow data from specimen
    ## and sampling data
    spec.data <- spec.data[spec.data$site %in% natural.sites,]
    spec.data <- spec.data[spec.data$site %in% names(d.area),]
    spec.data <- spec.data[spec.data$site %in% rownames(flower.mat),]
    spec.data$site <- as.character(spec.data$site)
    null.mat <- null.mat[rownames(null.mat) %in%
                         natural.sites,]
    null.mat <- null.mat[rownames(null.mat) %in%
                         names(d.area),]
    null.mat <- null.mat[rownames(null.mat) %in%
                         rownames(flower.mat),]

    ## make data array
    pollinator.id <- id(spec.data$pollinator)
    occ.mat <- make.mats(pollinator.id,
                         null.mat,
                         pollinator=as.vector(spec.data$pollinator),
                         var1=as.vector(spec.data$site),
                         var2=spec.data$date)
    sites <- rownames(occ.mat[[1]])
    dates <- colnames(occ.mat[[1]])
    species <- names(occ.mat)
    mat <- array(unlist(occ.mat),
                 dim=c(dim(occ.mat[[1]]),
                       length(occ.mat)))
    dimnames(mat) <- list(site=sites, date=dates, species=species)

    mm <- make.mat(mat, threshold, nzero)
    mat <- mm$mat

    ## make 4D matrix
    occ.mat <- lapply(1:dim(mat)[2], function(x) mat[,x,])
    occ.mat.split <- split(occ.mat, dimnames(mat)$date)
    yr.table <- table(dimnames(mat)$date)
    X <- array(NA, dim=c(dim(mat)[1], length(yr.table),
                         max(yr.table), dim(mat)[3]))
    dimnames(X) <- list(site=dimnames(mat)$site,
                        year=unique(dimnames(mat)$date),
                        rep=1:max(yr.table),
                        species=dimnames(mat)$species)

    null.mat <- matrix(NA, dim(mat)[1], dim(mat)[3],
                       dimnames=dimnames(X)[c('site', 'species')])
    f <- function(i) {
        missing <- max(yr.table) - yr.table[i]
        if(missing == 0) return(occ.mat.split[[i]])
        c(occ.mat.split[[i]], lapply(1:missing, function(x) null.mat))
    }
    tmp <- lapply(seq_along(yr.table), f)

    for(i in 1:length(yr.table))
        for(j in 1:max(yr.table))
            X[,i,j,] <- tmp[[i]][[j]]

    ## only keep sites with some positive number of reps
    no.reps <- apply(mat, 1, function(x) sum(x >= 0, na.rm=TRUE))==0
    X <- X[!no.reps,,,,drop=FALSE]

    ## create an array of day of the year
    make.date.mat <- function(dates) {
        date.mat <- array(NA, dim=dim(X)[1:3],
                          dimnames=dimnames(X)[1:3])
        for(i in seq_along(dates)) {
            year <- as.numeric(format(as.Date(dates[[i]],
                                              format='%Y-%m-%d'),
                                      format = '%Y'))
            lengths <- rle(as.vector(year))$lengths
            ind <- cbind(rep(i, sum(lengths)),
                         match(year, dimnames(X)$year),
                         as.vector(unlist(sapply(lengths,
                                                 seq_len))))
            date.mat[ind] <- strptime(dates[[i]],
                                      '%Y-%m-%d')$yday+1
        }
        date.mat
    }
    date.occ.mat <- lapply(mm$dates, make.date.mat)

    date.mat <- array(unlist(date.occ.mat),
                      dim=c(dim(date.occ.mat[[1]]), length(date.occ.mat)))
    dimnames(date.mat) <- dimnames(X)

    ## function to re-arrange replicate dimension
    compress <- function(x) {
        if(!any(is.na(x))) return(x)
        return(c(x[!is.na(x)], x[is.na(x)]))
    }
    X <- aperm(apply(X, c(1,2,4), compress), c(2,3,1,4))
    names(dimnames(X)) <- c("site", "year", "rep", "species")
    date.mat <- aperm(apply(date.mat, c(1,2,4), compress), c(2,3,1,4))

    ## extract traits, scale
    these.traits <- traits[, c(col.name.trait1, col.name.trait2) ]
    rownames(these.traits) <- traits$GenusSpecies
    these.traits <- these.traits[rownames(these.traits) %in%
                                 dimnames(X)$sp,]
    ## log and standardize trait data
    these.traits <- apply(these.traits, 2,
                          function(x) standardize(log(x + 1)))
    any.na.trait <- apply(these.traits, 1,
                          function(x) !any(is.na(x)))
    traits1 <- these.traits[any.na.trait,1]
    traits2 <- these.traits[any.na.trait,2]
    X <- X[,,,any.na.trait,drop=FALSE]
    date.mat <- date.mat[,,,any.na.trait,drop=FALSE]

    ## drop last year of flower data
    flower.mat <- flower.mat[rownames(flower.mat) %in%
                             dimnames(X)$site,]
    flower.mat <- flower.mat[, as.numeric(colnames(flower.mat)) <
                               max(dimnames(X)$year)]
    flower.mat <- standardize(flower.mat)

    ## drop last year of natural data, log, standardize
    natural.mat <- natural.mat[names(natural.mat) %in%
                               dimnames(X)$site]
    natural.mat <- log(natural.mat + 1)
    natural.mat <- standardize(natural.mat)

    ## specify the initial values
    ## if a species was present, NA, if NA or 0, 1
    zinits <- apply(X, c(1,2,4),
                    function(x) (sum(x,na.rm=TRUE) > 0)*1)
    Z <- zinits
    ## data, all 0s should be NAs, i.e., should be sampled
    Z[which(zinits != 1)] <- NA

    ## don't sample data we know...
    zinits[zinits == 1] <- NA

    nsp <- dim(X)[4]
    inits <- c(list(Z=zinits),
               getInits(nsp))
    constants <-  list(nrep=apply(X, c(1,2,4),
                                  function(x) sum(x >= 0,na.rm=TRUE)),
                       nsp=dim(X)[4],
                       nsite=dim(X)[1],
                       nyear=dim(X)[2],
                       max.nreps = dim(X)[3])

    ## create polynomial terms
    day <- day.2 <- date.mat
    poly.terms <- poly(date.mat[!is.na(date.mat)], 2)
    day[!is.na(day)]  <- standardize(poly.terms[,1])
    day.2[!is.na(day.2)]  <- standardize(poly.terms[,2])

    data <- list(Z=Z,
                 X=X,
                 k=traits1,
                 B=traits2,
                 day=day,
                 day.2=day.2,
                 HRarea=d.area[names(d.area) %in% dimnames(X)$site][dimnames(X)$site],
                 natural=natural.mat[dimnames(X)$site],
                 fra = flower.mat[dimnames(X)$site,])

    monitors <- getParams()

    my.inits <- inits
    model.input <- list(data=data,
                        monitors=monitors,
                        constants=constants,
                        inits=my.inits)
    return(model.input)
}

samp2site.ypr <- function(site, yr, abund) {
    x <- tapply(abund, list(site=site,yr=yr), sum, na.rm=TRUE)
    return(x)
}


prepModel <- function(){
    if(filtering){
        source('src/dynamicOcc.R')
        model.input$data$Z <- NULL
        model.input$inits$Z <- NULL
        ## We do not want any X element equal to NA or they will not be
        ## considered data and will be sampled.
        model.input$data$X[ is.na(model.input$data$X) ] <- -1000
    }
    testOccData()
    return(model.input)
}



## ## return the closest syd in a data-set to argument syd
## get.close.syd <- function(syd, d, threshold=10) {
##     d.sub <- d[d$Site==get.site(syd) & d$Year ==as.character(get.year(syd)),]
##     days.diff <- abs(as.numeric(d.sub$Day)-get.day(syd))

##     if(length(days.diff)==0) {
##         cat(sprintf('no match for %s\n', syd))
##         return(NA)
##     }
##     if(min(days.diff)>threshold) {
##         cat(sprintf('date discrepancy > %s days: ', threshold))
##         cat(sprintf('%s\n', syd))
##     }
##     d.sub$syd[which.min(days.diff)]
## }

## get.site <- function(syd)
##     as.vector(sapply(syd, function(x)
##         strsplit(x, ';')[[1]][1]))
## get.year <- function(syd)
##     as.vector(sapply(syd, function(x)
##         as.numeric(strsplit(x, ';')[[1]][2])))
## get.day <- function(syd)
##     as.vector(sapply(syd, function(x)
##         as.numeric(strsplit(x, ';')[[1]][3])))



