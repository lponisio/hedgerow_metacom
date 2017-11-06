
makeNetworkFig <- function(spec,
                           sys=NULL,
                           lat.long=NULL,
                           rows,
                           cols,
                           nets,
                           weight.edge.colors=FALSE,
                           cols.vertex,
                           cols.edges,
                           vertex.names = c("control", "mature"),
                           natural.cover=NULL,
                           site.ave=NULL, ...){
    expandNets <- function(sub.mat, all.mat){
        sub.mat <- sub.mat[rownames(sub.mat) %in% rownames(all.mat),]
        sub.mat <- sub.mat[,colnames(sub.mat) %in% colnames(all.mat)]
        indx1 <- match(rownames(sub.mat), rownames(all.mat))
        indx2 <-  match(colnames(sub.mat), colnames(all.mat))
        all.mat[indx1,
                indx2] <- sub.mat
        return(all.mat)
    }
    ## create an expanded matrix of all the pollinators/sites
    pols <- unique(spec[, cols])
    sites <- unique(spec[, rows])

    all.pp <- matrix(NA, nrow=length(sites),
                     ncol=length(pols))
    rownames(all.pp) <- sites
    colnames(all.pp) <- pols

    nets.all <- simplify2array(rapply(nets, expandNets,
                                      all.mat=all.pp,
                                      how="replace"))

    mean.net.year <- apply(nets.all, c(1,2), sum, na.rm=TRUE)
    mean.net.year[is.na(mean.net.year)] <- 0
    mean.net.year[mean.net.year > 0] <- 1
    dist.year <- designdist(mean.net.year, method="J")
    dist.year <- as.matrix(dist.year)
    if(!is.null(site.ave)){
        importance <- site.ave
    } else{
        importance <-  (rowSums(dist.year)/
                        max(rowSums(dist.year)))*2
    }
    ## zero out duplicate edges
    dist.year[lower.tri(dist.year, diag=FALSE)] <- 0

    ## great graph
    gs <- graph.incidence(dist.year, weighted=TRUE)

    ## color vertices by site stype, weight by the number of shared
    ## species,  add lat long positions
    V(gs)$size <- importance
    E(gs)$width <- (E(gs)$weight/max(E(gs)$weight))*1.5

    site.status <- spec$SiteStatus[match(V(gs)$name,
                                         spec$Site)]
    v.boarders <- rep("black", length(importance))

    if(!is.null(lat.long)){
        site.status <- spec$SiteStatus[match(V(gs)$name,
                                             spec$Site)]
        names(site.status) <- V(gs)$name
        names(cols.vertex) <- vertex.names
        V(gs)$color <- cols.vertex[site.status]
        poss <- match(V(gs)$name, rownames(lat.long))
        gs <- delete_vertices(gs, V(gs)[is.na(poss)])
        poss <- poss[!is.na(poss)]
        V(gs)$long <- lat.long[,1][poss]
        V(gs)$lat <- lat.long[,2][poss]
        latex <- matrix(c(V(gs)$long, V(gs)$lat), ncol=2)
        colnames(latex) <- c("long", "lat")
        gs$layout <- latex

        edge.names <-  get.edgelist(gs)
        status.edges <- cbind(site.status[edge.names[,1]],
                              site.status[edge.names[,2]])
        edge.types <- apply(status.edges, 1, function(x) paste(x,
                                                               collapse="_"))
        names(cols.edges) <- unique(edge.types)
        E(gs)$color <- cols.edges[edge.types]

        dims <- bbox(sys)
        plot(sys)
        rect(xleft=dims[1,1],ybottom=dims[2,1],
             xright=dims[1,2],ytop=dims[2,2],
             col= rgb(1,1,1, alpha=0.3))
        if(!is.null(natural.cover)){
            natural.cover <- crop(natural.cover,  dims)
            poly.cols <- add.alpha("black",
                                   alpha=0.2)
            plot(natural.cover,
                 col=poly.cols, border=poly.cols, add=TRUE)
        }

    }


    if(weight.edge.colors){
        ## Color scaling function
        c_scale <- colorRamp(c('white','black'))

        ## Applying the color scale to edge weights.
        ## rgb method is to convert colors to a character vector.
        E(gs)$color <- apply(c_scale(E(gs)$weight/max(E(gs)$weight)), 1,
                             function(x){
                                 rgb(x[1]/255,x[2]/255,x[3]/255,
                                     alpha=0.7)
                             })
    }


    plot(gs, vertex.label="",
         vertex.frame.color=v.boarders,
         edge.curved=0.4,...)
    legend("bottomright", legend=c("Surveyed field margin", "Surveyed hedgerow",
                                   "Unsurveyed hedgerow"),
           col=c(cols.vertex, "goldenrod2"),
           pch=c(16,16,15), cex=0.75, bg="white", inset=c(0.036,0.04))
    legend("bottomright", legend=c("Surveyed field margin", "Surveyed hedgerow",
                                   "Unsurveyed hedgerow"),
           col="black", pch=c(1,1,0), cex=0.75,inset=c(0.036,0.04))

    if(!is.null(lat.long)){
        points(latex, col=V(gs)$color, pch=16, cex=importance)
        points(latex, col="black", pch=1, cex=importance)
        new.sites <-
            rownames(lat.long)[!rownames(lat.long) %in%
                               unique(spec$Site)]
        points(lat.long[new.sites,], col="black", pch=0,
               cex=0.3)
        points(lat.long[new.sites,], col="goldenrod2", pch=15,
               cex=0.25)

    }
}

plotAllStatuses <- function(){
    cols.vertex <- brewer.pal(11, 'RdGy')[c(1,10)]
    cols.edges <- add.alpha(brewer.pal(11, 'RdGy')[c(6,10,1,6)], alpha=0.5)
    makeNetworkFig(spec, sys, lat.long=lat.long,
                   rows="Site",
                   cols="GenusSpecies",
                   nets=nets.year,
                   natural.cover=landcover.nat,
                   add=TRUE,
                   rescale = FALSE,
                   site.ave=this.site.ave,
                   cols.vertex=cols.vertex,
                   cols.edges=cols.edges)

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

