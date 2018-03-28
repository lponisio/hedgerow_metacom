## ************************************************************
## prep and format data
## ************************************************************
samp2site.ypr <- function(site, yr, abund) {
    x <- tapply(abund, list(site=site,yr=yr), sum, na.rm=TRUE)
    return(x)
}


## return the closest syd in a data-set to argument syd
get.close.syd <- function(syd, d, threshold=10) {
    d.sub <- d[d$Site==get.site(syd) & d$Year ==as.character(get.year(syd)),]
    days.diff <- abs(as.numeric(d.sub$Day)-get.day(syd))

    if(length(days.diff)==0) {
        cat(sprintf('no match for %s\n', syd))
        return(NA)
    }
    if(min(days.diff)>threshold) {
        cat(sprintf('date discrepancy > %s days: ', threshold))
        cat(sprintf('%s\n', syd))
    }
    d.sub$syd[which.min(days.diff)]
}

get.site <- function(syd)
    as.vector(sapply(syd, function(x)
        strsplit(x, ';')[[1]][1]))
get.year <- function(syd)
    as.vector(sapply(syd, function(x)
        as.numeric(strsplit(x, ';')[[1]][2])))
get.day <- function(syd)
    as.vector(sapply(syd, function(x)
        as.numeric(strsplit(x, ';')[[1]][3])))


prepOccModelInput <- function(nzero, ## if augmenting data
                              threshold, ## min times a pollinator is seen
                              save.dir, ## directory to save data
                              spec, ## specimen data
                              sr.sched, ## sampling schedule
                              traits, ## dataframe with traits
                              col.name.trait1, ## colname of one trait
                              col.name.trait2, ## colname of second trait
                              HRarea, ## vector of hedgerow area prox
                              buffer=NULL, ## if using buffers, what
                              ## is the d of the buffer of interest
                              natural.mat=NULL, ## site by year
                              ## natural area proximity
                              natural.decay,## site by decay matrix natural
                              veg, ## floral availability matrix
                              col.name.div.type="Div", ## colname of
                              ## the metric of floral availability
                              load.inits.name='nimble', ## file name
                              ## of the model to load inits from
                              jags=FALSE, ## prep input for jags?
                              model.type="allInt",
                              kremen.digitize=FALSE,
                              HR.decay=TRUE,
                              data.subset="all",
                              raw.flower=FALSE) {

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
    if(!is.null(buffer)){
        d.area <- HRarea[, buffer]
    } else if(HR.decay){
        hr.mat <- HRarea[, natural.decay]
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

    if(!raw.flower){
        flower.mat <- samp2site.ypr(site=veg$Site,
                                    yr=veg$Year,
                                    abund=veg[, col.name.div.type])
    } else{
        sr.sched <- sr.sched[sr.sched$NetPan == "net",]
        sr.sched$Day <- strftime(sr.sched$Date, format="%j")
        sr.sched$Year <- format(as.Date(sr.sched$Date),"%Y")
        veg$Day <- strftime(veg$Date, format="%j")
        veg$Year <- format(as.Date(veg$Date),"%Y")
        veg$syd <- paste(veg$Site, veg$Year, veg$Day,
                         sep=";")
        syd <- paste(sr.sched$Site, sr.sched$Year, sr.sched$Day,
                     sep=";")
        close.veg.data <- sapply(syd, get.close.syd, d=veg,
                                 threshold=30)
        ## WIP


    }
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

    ## drop last year of flower data,  don't standardize, dont log!
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

    if(model.type == "no_noncrop"){
        data$natural <- NULL
        inits <- inits[!grepl("nat", names(inits))]
        monitors <- monitors[!grepl("nat", monitors)]
    }

    if(jags){
        my.inits <- function() inits
    } else {
        my.inits <- inits
    }
    save.path <- file.path(save.dir,
                           sprintf('%s-%s-%s-%s.RData',
                                   data.subset, threshold,
                                   nzero, natural.decay))

    model.input <- list(data=data,
                        monitors=monitors,
                        constants=constants,
                        inits=my.inits)
    if(!dir.exists(save.dir)) {
        cat(paste("Needed dir", save.dir, "does not exist. OK to create? (Type 'yes' if ok.)"))
        okToMakeDir <- readlines()
        if(!identical(okToMakeDir, "yes"))
            stop("Stopping because permission to make save.dir denied.")
        dir.create(save.dir, showWarnings = FALSE)
    }
    save(model.input, file=save.path)
    return(model.input)
}

## ************************************************************
## specify the parameters to be monitored
## ************************************************************

getParams <- function(){
    c('mu.p.0',
      'sigma.p.0',
      'mu.p.day.1',
      'sigma.p.day.1',
      'mu.p.day.2',
      'sigma.p.day.2',

      'mu.phi.0',
      'sigma.phi.0',
      'mu.phi.hr.area',
      'sigma.phi.hr.area',
      'mu.phi.nat.area',
      'sigma.phi.nat.area',
      'mu.phi.fra',
      'sigma.phi.fra',
      'phi.k',

      'phi.B',
      'phi.nat.area.fra',
      'phi.hr.area.fra',
      'phi.nat.area.k',
      'phi.hr.area.k',
      'phi.nat.area.B',
      'phi.hr.area.B',

      'mu.gam.0',
      'sigma.gam.0',
      'mu.gam.hr.area',
      'sigma.gam.hr.area' ,
      'mu.gam.nat.area',
      'sigma.gam.nat.area',
      'mu.gam.fra',
      'sigma.gam.fra',

      'gam.k',
      'gam.B',
      'gam.hr.area.fra',
      'gam.nat.area.fra',
      'gam.hr.area.k',
      'gam.nat.area.k',
      'gam.hr.area.B',
      'gam.nat.area.B'


      ## ## site level effects
      ## 'phi.nat.area',
      ## 'phi.hr.area',
      ## 'phi.hr.area.fra',
      ## 'phi.nat.area.fra',
      ## 'phi.hr.area.k'

      )
}

getInits <- function(nsp){
    list(mu.p.0 = rnorm(1),
         sigma.p.0 = runif(1),
         mu.p.day.1 = rnorm(1),
         sigma.p.day.1 = runif(1),
         mu.p.day.2 = rnorm(1),
         sigma.p.day.2= runif(1),

         mu.phi.0 = rnorm(1),
         sigma.phi.0 = runif(1),

         mu.phi.hr.area = rnorm(1),
         sigma.phi.hr.area = runif(1),
         mu.phi.nat.area = rnorm(1),
         sigma.phi.nat.area = runif(1),
         mu.phi.fra = rnorm(1),
         sigma.phi.fra = runif(1),

         mu.gam.0 = rnorm(1),
         sigma.gam.0 = runif(1),

         mu.gam.hr.area = rnorm(1),
         sigma.gam.hr.area = runif(1),
         mu.gam.nat.area = rnorm(1),
         sigma.gam.nat.area = runif(1),
         mu.gam.fra = rnorm(1),
         sigma.gam.fra = runif(1),

         p.0 = rnorm(nsp),
         p.day.1 = rnorm(nsp),
         p.day.2 = rnorm(nsp),

         phi.0 = rnorm(nsp),
         phi.hr.area = rnorm(nsp),
         phi.nat.area = rnorm(nsp),
         phi.fra = rnorm(nsp),

         phi.k = rnorm(1),
         phi.B = rnorm(1),
         phi.hr.area.fra = rnorm(1),
         phi.nat.area.fra = rnorm(1),
         phi.nat.area.k = rnorm(1),
         phi.hr.area.k = rnorm(1),
         phi.nat.area.B = rnorm(1),
         phi.hr.area.B = rnorm(1),

         gam.0 = rnorm(nsp),
         gam.hr.area = rnorm(nsp),
         gam.nat.area = rnorm(nsp),
         gam.fra = rnorm(nsp),

         gam.k = rnorm(1),
         gam.B = rnorm(1),
         gam.hr.area.fra = rnorm(1),
         gam.nat.area.fra = rnorm(1),
         gam.hr.area.k = rnorm(1),
         gam.nat.area.k = rnorm(1),
         gam.hr.area.B = rnorm(1),
         gam.nat.area.B = rnorm(1)
         )

}
