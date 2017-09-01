## ************************************************************
## prep and format data
## ************************************************************
samp2site.ypr <- function(site, yr, abund) {
    x <- tapply(abund, list(site=site,yr=yr), sum, na.rm=TRUE)
    return(x)
}

prepOccModelInput <- function(nzero, ## if augmenting data
                              threshold, ## min times a pollinatoris seen
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
                              natural.decay, ## site by decay matrix
                              veg, ## floral availability matrix
                              w.ypr=FALSE, ## include ypr in model?
                              load.inits=TRUE, ## load inits from a
                              ## prior model?
                              init.name="nimble", ## name of the
                              ## proior model to load inits from
                              col.name.div.type="Div", ## colname of
                              ## the metric of floral availability
                              load.inits.name='nimble', ## file name
                              ## of the model to load inits from
                              jags=FALSE) { ## prep input for jags?

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

    ## natural cover, scale by year
    ## for the kremen digitized data, the antural data is year
    ## specific so a matrix, otherwise it is a vector
    if(is.null(natural.mat)){ ## for kremen digitized data
        natural.mat <- samp2site.spp(natural$Site,
                                     natural$Year, natural[,
                                                           natural.decay])
        natural.mat <- apply(natural.mat, 2, standardize)
        natural.sites <- rownames(natural.mat)
    } else{ ## for yolo landcover data
        natural.mat <- natural.mat[, natural.decay]
        natural.mat <- natural.mat[!is.na(natural.mat)]
        natural.mat <- standardize(natural.mat)
        natural.sites <- names(natural.mat)
    }

    ## area of hedgerow, scaled
    if(!is.null(buffer)){
        d.area <- HRarea[, buffer]

    } else{
        d.area <- HRarea
        buffer <- "all"
    }
    d.area <- standardize(d.area[names(d.area) %in%
                                 natural.sites])

    ## veg diversity
    veg$Site <- gsub(":.*", "", veg$Site)
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
    no.reps <- apply(mat, 1, function(x) sum( x>= 0, na.rm=TRUE))==0
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
    these.traits <- apply(these.traits, 2, standardize)
    any.na.trait <- apply(these.traits, 1, function(x) !any(is.na(x)))
    traits1 <- these.traits[any.na.trait,1]
    traits2 <- these.traits[any.na.trait,2]
    X <- X[,,,any.na.trait,drop=FALSE]
    date.mat <- date.mat[,,,any.na.trait,drop=FALSE]

    ## drop last year of flower data, standardize
    flower.mat <- flower.mat[rownames(flower.mat) %in%
                             dimnames(X)$site,]
    flower.mat <- flower.mat[, colnames(flower.mat) !=
                               max(dimnames(X)$year)]
    flower.mat <- apply(flower.mat, 2, standardize)


    ## drop last year of natural data
    if(is.null(natural.mat)){
        natural.mat <- natural.mat[, colnames(natural.mat) !=
                                     max(dimnames(X)$year)]
        natural.mat <- natural.mat[rownames(natural.mat) %in%
                                   dimnames(X)[[1]],]
    } else {
        natural.mat <- natural.mat[names(natural.mat) %in%
                                   dimnames(X)[[1]]]
    }
    site.status <- unique(cbind(Site=spec$Site,
                                Year=spec$Year,
                                ypr=spec$ypr,
                                SiteStatus=spec$SiteStatus))

    ## convert occurrence data into a site by species matrix
    site.status.mat <- samp2site.ypr(site=site.status[, "Site"],
                                     yr=site.status[, "Year"],
                                     abund=as.numeric(site.status[,
                                                                  "ypr"]))

    mature <- site.status[,"Site"][site.status[, "SiteStatus"] ==
                                   "mature"]
    maturing <- site.status[,"Site"][site.status[, "SiteStatus"] ==
                                     "maturing"]
    site.status.mat[rownames(site.status.mat) %in% mature,] <-
        rep(10, ncol(site.status.mat))

    controls <- apply(site.status.mat, 1, sum, na.rm=TRUE)
    controls <- names(controls)[controls == 0]

    site.status.mat[rownames(site.status.mat) %in% controls,] <-
        rep(0, ncol(site.status.mat))

    for(site in maturing){
        this.site <- site.status.mat[site, ]
        ypr1 <- names(this.site)[this.site == min(this.site, na.rm=TRUE) &
                                 !is.na(this.site)]
        this.site <-
            as.numeric(names(this.site)) -
            (as.numeric(ypr1) -1)
        this.site[this.site < 0] <- 0
        site.status.mat[site, ] <- this.site
    }

    site.status.mat <- site.status.mat[rownames(site.status.mat) %in%
                                       dimnames(X)$site,
                                       colnames(site.status.mat) !=
                                       max(dimnames(X)$year)]

    ## specify the initial values
    ## if a species was present, NA, if NA or 0, 1
    zinits <- apply(X, c(1,2,4),
                    function(x) (sum(x,na.rm=TRUE) > 0)*1)
    ## zinits[apply(X, c(1,2,4), function(x) !any(!is.na(x)))] <- NA
    Z <- zinits
    Z[which(zinits != 1)] <- NA
    zinits[Z == 1] <- NA

    nsp <- dim(X)[4]
    inits <- c(list(Z=zinits),
                  getInits(nsp))

    ## load initial conditions from a previous run
    if(load.inits){
        dd <- load(file=file.path(save.dir, sprintf('runs/%s.Rdata',
                                                    load.inits.name)))
        load(file=file.path(save.dir, sprintf('runs/%s.Rdata',
                                              load.inits.name)))

        dd <- eval(parse(text =dd))

        inits <- c(dd[[1]]$summary[init.name,
                                   "median",], inits)
        inits <- inits[!duplicated(names(inits))]
    }

    save.path <- file.path(save.dir,
                           sprintf('%s-%s-%s.RData',
                                   threshold, nzero, buffer))
    constants <-  list(nrep=apply(X, c(1,2,4),
                                  function(x) sum(x >= 0,na.rm=TRUE)),
                       nsp=dim(X)[4],
                       nsite=dim(X)[1],
                       nyear=dim(X)[2])

    ## create polynomial terms
    day <- day.2 <- date.mat
    poly.terms <- poly(date.mat[!is.na(date.mat)], 2)
    day[!is.na(day)]  <- standardize(poly.terms[,1])
    day.2[!is.na(day.2)]  <- standardize(poly.terms[,2])

    data <- list(Z=Z,
                 X=X,
                 traits1=traits1,
                 traits2=traits2,
                 day=day,
                 day.2=day.2,
                 ypr = site.status.mat,
                 HRarea=d.area[names(d.area) %in% dimnames(X)[[1]]],
                 natural=natural.mat,
                 fra = flower.mat)

    ## remove ypr and inits if not in model
    if(!w.ypr){
        data$ypr <- NULL
        inits <- inits[!grepl("ypr", names(inits))]
    }

    if(jags){
        my.inits <- function() inits
    } else {
        my.inits <- inits
    }

    model.input <- list(data=data,
                        monitors=getParams(),
                        constants=constants,
                        inits=my.inits)
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
      'mu.phi.traits1',
      'sigma.phi.traits1',
      'mu.phi.traits2',
      'sigma.phi.traits2',
      'mu.phi.nat.area.fra',
      'sigma.phi.nat.area.fra',
      'mu.phi.hr.area.fra',
      'sigma.phi.hr.area.fra',
      'mu.phi.nat.area.traits1',
      'sigma.phi.nat.area.traits1',
      'mu.phi.hr.area.traits1',
      'sigma.phi.hr.area.traits1',


      'mu.gam.0',
      'sigma.gam.0',
      'mu.gam.hr.area',
      'sigma.gam.hr.area',
      'mu.gam.nat.area',
      'sigma.gam.nat.area',
      'mu.gam.fra',
      'sigma.gam.fra',
      'mu.gam.traits1',
      'sigma.gam.traits1',
      'mu.gam.traits2',
      'sigma.gam.traits2',
      'mu.gam.traits1.fra',
      'sigma.gam.traits1.fra',
      'mu.phi.traits1.fra',
      'sigma.phi.traits1.fra',
      'mu.gam.hr.area.fra',
      'sigma.gam.hr.area.fra',
      'mu.gam.nat.area.fra',
      'sigma.gam.nat.area.fra',
      'mu.gam.hr.area.traits1',
      'sigma.gam.hr.area.traits1',
      'mu.gam.nat.area.traits1',
      'sigma.gam.nat.area.traits1',

      'phi.site.mean',
      'gam.site.mean',
      'phi.sp.mean',
      'gam.sp.mean')
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
      mu.phi.traits1 = rnorm(1),
      sigma.phi.traits1 = runif(1),
      mu.phi.traits2 = rnorm(1),
      sigma.phi.traits2 = runif(1),
      mu.phi.nat.area.fra = rnorm(1),
      sigma.phi.nat.area.fra = runif(1),
      mu.phi.hr.area.fra = rnorm(1),
      sigma.phi.hr.area.fra = runif(1),
      mu.phi.nat.area.traits1 = rnorm(1),
      sigma.phi.nat.area.traits1 = runif(1),
      mu.phi.hr.area.traits1 = rnorm(1),
      sigma.phi.hr.area.traits1 = runif(1),

      mu.gam.0 = rnorm(1),
      sigma.gam.0 = runif(1),
      mu.gam.hr.area = rnorm(1),
      sigma.gam.hr.area = runif(1),
      mu.gam.nat.area = rnorm(1),
      sigma.gam.nat.area = runif(1),
      mu.gam.fra = rnorm(1),
      sigma.gam.fra = runif(1),
      mu.gam.traits1 = rnorm(1),
      sigma.gam.traits1 = runif(1),
      mu.gam.traits2 = rnorm(1),
      sigma.gam.traits2 = runif(1),
      mu.gam.traits1.fra = rnorm(1),
      sigma.gam.traits1.fra = runif(1),
      mu.phi.traits1.fra = rnorm(1),
      sigma.phi.traits1.fra = runif(1),
      mu.gam.hr.area.fra = rnorm(1),
      sigma.gam.hr.area.fra = runif(1),
      mu.gam.nat.area.fra = rnorm(1),
      sigma.gam.nat.area.fra = runif(1),
      mu.gam.hr.area.traits1 = rnorm(1),
      sigma.gam.hr.area.traits1 = runif(1),
      mu.gam.nat.area.traits1 = rnorm(1),
      sigma.gam.nat.area.traits1 = runif(1),

      p.0 = rnorm(nsp),
      p.day.1 = rnorm(nsp),
      p.day.2 = rnorm(nsp),

      phi.0 = rnorm(nsp),
      phi.hr.area = rnorm(nsp),
      phi.nat.area = rnorm(nsp),
      phi.fra = rnorm(nsp),
      phi.traits1 = rnorm(nsp),
      phi.traits2 = rnorm(nsp),
      phi.nat.area.fra = rnorm(nsp),
      phi.hr.area.fra = rnorm(nsp),
      phi.nat.area.traits1 = rnorm(nsp),
      phi.hr.area.traits1 = rnorm(nsp),

      gam.0 = rnorm(nsp),
      gam.hr.area = rnorm(nsp),
      gam.nat.area = rnorm(nsp),
      gam.fra = rnorm(nsp),
      gam.traits1 = rnorm(nsp),
      gam.traits2 = rnorm(nsp),
      gam.traits1.fra = rnorm(nsp),
      phi.traits1.fra = rnorm(nsp),
      gam.hr.area.fra = rnorm(nsp),
      gam.nat.area.fra = rnorm(nsp),
      gam.hr.area.traits1 = rnorm(nsp),
      gam.nat.area.traits1 = rnorm(nsp))
}
