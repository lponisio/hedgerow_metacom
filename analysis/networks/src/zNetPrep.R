
species.lev <- function(x){
    present.x <- bipartite::empty(x)
    sl <- specieslevel(present.x,
                       index=c("degree", "betweenness", "closeness"))
    if(ncol(x) != ncol(present.x)){
        not.present <- colnames(x)[!colnames(x) %in%
                                   rownames(sl$"higher level")]
        all.na.mat <- matrix(NA, ncol=ncol(sl$'higher level'),
                             nrow=length(not.present))
        rownames(all.na.mat) <- not.present
        colnames(all.na.mat) <- colnames(sl$'higher level')
        sl$'higher level' <- rbind(sl$'higher level', all.na.mat)
        sl$'higher level' <-
            sl$'higher level'[order(rownames(sl$'higher level')),]
    }

    if(nrow(x) != nrow(present.x)){
        not.present <- rownames(x)[!rownames(x) %in%
                                   rownames(sl$"lower level")]
        all.na.mat <- matrix(NA, ncol=ncol(sl$'lower level'),
                             nrow=length(not.present))
        rownames(all.na.mat) <- not.present
        colnames(all.na.mat) <- colnames(sl$'lower level')
        sl$'lower level' <- rbind(sl$'lower level', all.na.mat)
        sl$'lower level' <-
            sl$'lower level'[order(rownames(sl$'lower level')),]
    }
    return(sl)
}

calcMetrics <- function(j, z.array, by.year){
    if(by.year){
        this.year <- z.array[,j,]
    } else{
        this.year <- z.array[j,,]
    }
    specs.year <- species.lev(this.year)
    ll <- specs.year$`lower level`
    hl <- specs.year$`higher level`
    return(list(ll=ll, hl=hl))
}


calcYearNet <- function(i, these.samples, model.input){
    print(i)
    z.array <- array(these.samples[i,],
                     dim=dim(model.input$data$Z),
                     dimnames=dimnames(model.input$data$Z))
    metrics.yr <- lapply(1:dim(z.array)[2], calcMetrics, z.array,
                           by.year=TRUE)
    metrics.site <- lapply(1:dim(z.array)[1], calcMetrics, z.array,
                             by.year=FALSE)
    return(list(site= lapply(metrics.yr, function(x) x[["ll"]]),
                sp=lapply(metrics.yr, function(x) x[["hl"]]),
                year=lapply(metrics.site, function(x) x[["hl"]])))
}


getPosteriorNet <- function(j, nets, model.input, site.sp){
    print(j)
    this.year <- lapply(nets, function(x) as.matrix(x[[j]]))
    this.year.array <- simplify2array(this.year)

    year.mean <- apply(this.year.array, c(1,2), mean,
                       na.rm=TRUE)
    year.sd <- apply(this.year.array, c(1,2), sd,
                     na.rm=TRUE)
    colnames(year.sd) <- paste0("sd.", colnames(year.sd))
    year <- as.data.frame(cbind(year.mean, year.sd))
    if(site.sp == "site"){
        year$Year <- dimnames(model.input$data$Z)[[2]][j]
        year$Site <- rownames(year)
    } else if(site.sp == "sp"){
        year$Year <- dimnames(model.input$data$Z)[[2]][j]
        year$GenusSpecies <- rownames(year)
    } else if(site.sp == "year"){
        year$Site <- dimnames(model.input$data$Z)[[1]][j]
        year$GenusSpecies <- rownames(year)
    }
    rownames(year) <- NULL
    return(year)
}
