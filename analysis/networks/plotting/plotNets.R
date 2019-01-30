## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')

## source('src/initialize.R')
source('src/misc.R')
library(bipartite)
library(igraph)
f.path <- '../../data/networks'
load(file=file.path(f.path, 'all_networks_years.Rdata'))
load(file=file.path(f.path, 'years_networks.Rdata'))
load(file=file.path(f.path, 'sites_networks.Rdata'))

## nets.site is within a site, across years
## nets.year is within a year, across sites

fig.path <- "../../../hedgerow_metacom_saved/networks/figures"

plotNet <- function(nets, type, index){
    par(mar=c(0,0.2,1,0.5), oma=c(0,0,0,0))
    rows <- ifelse((length(nets) %% 2) == 0, 2, 3)
    layout(matrix(1:length(nets), nrow=rows, byrow=TRUE))
    for(j in 1:length(nets)){
        g <- nets[[j]]
        gs <- graph.incidence(g, weighted=TRUE)
        cols <- c(rep("darkolivegreen", length(rownames(g))),
                  rep("gold", length(colnames(g))))
        if(index=="degree"){
            importance <-  (c(rowSums(g) +0.1, colSums(g) +
                                               0.1)/sum(g))
        } else{
            metrics <- specieslevel(g, index=index)
            metrics <- sapply(metrics, function(x){
                x[,paste0("weighted.",index)]
            })
            importance <- c(metrics$'lower level',
                            metrics$'higher level')
            importance <- importance/max(importance)

        }
        V(gs)$color <- cols
        v.labs <- names(V(gs))
        v.labs[nchar(v.labs) > 9] <- ""
        V(gs)$size <- importance*50
        v.boarders <- rep("black", length(importance))
        if(type=="year"){
            V(gs)$size[v.labs != ""] <- 0
            V(gs)$color[v.labs != ""] <- "white"
            v.boarders[v.labs != ""] <- "white"
            dists <- 2
        } else{
            v.labs <- ""
            dists <- 0
        }
        E(gs)$color <- add.alpha("grey")
        ## E(gs)$width <- (E(gs)$weight/sum(E(gs)$weight))*10
        gs$layout <- layout_in_circle

        plot.igraph(gs, vertex.frame.color=v.boarders,
                    vertex.label.cex=0.5,
                    vertex.label=v.labs,
                    vertex.label.dist=dists)
        mtext(names(nets)[j], 3, cex=0.5)
    }
}



plotSiteNets <- function(){plotNet(nets.site, "year", "degree")}

pdf.f(plotSiteNets, file=file.path(fig.path,
                                   sprintf("siteNets%s.pdf",
                                           "degree")),
      width=8, height=4)

plotSiteNets <- function(){plotNet(nets.site, "year", "betweenness")}

pdf.f(plotSiteNets, file=file.path(fig.path,
                                   sprintf("siteNets%s.pdf",
                                           "betweenness")),
      width=8, height=4)



plotYearNets <- function(){plotNet(nets.year, "site", "degree")}
pdf.f(plotYearNets, file=file.path(fig.path,
                                   sprintf("yearNets%s.pdf",
                                           "degree")),
      width=7, height=3)

plotYearNets <- function(){plotNet(nets.year, "site", "betweenness")}
pdf.f(plotYearNets, file=file.path(fig.path,
                                   sprintf("yearNets%s.pdf",
                                           "betweenness")),
      width=7, height=3)
