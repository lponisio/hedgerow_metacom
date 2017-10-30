
makePredictData <- function(dd, y, mods, time=TRUE){
    dd <- cbind(dd, 0)
    colnames(dd) <- c(colnames(dd)[colnames(dd) != "0"], y)
    dd.pi.space <- predict.int(mod= mods[["space"]][[y]],
                               dd=dd,
                               y=y,
                               family="gaussian")
    if(time){
        dd.pi.time <- predict.int(mod= mods[["time"]][[y]],
                                  dd=dd,
                                  y=y,
                                  family="gaussian")
        return(list(time=dd.pi.time,
                    space=dd.pi.space))
    }else{
        return(list(space=dd.pi.space))
    }

}

plot.panels <- function(){
    f <- function(){
        col.lines <- "black"
        col.fill <- add.alpha(col.lines, alpha=0.2)
        col.white <- add.alpha("white", alpha=0)
        names(col.lines) <- names(col.fill) <- "all"
        layout(matrix(1:4, nrow=2, byrow=TRUE))
        par(oma=c(4, 7, 3, 1),
            mar=c(1, 0, 0, 1.5), cex.axis=1.2)

        ys <- c("k", "betweenness")
        ylabs <- c("Degree (z)", "Centrality")

        dd.pi.degree.degree <- makePredictData(r.degree.dd, ys[1], pol.mods)
        dd.pi.degree.between <- makePredictData(r.degree.dd, ys[2], pol.mods)

        ## space degree-degree
        plot.panel(new.dd=dd.pi.degree.degree$space,
                   dats=pol.specs[["space"]],
                   y1=ys[1],
                   y2= range(c(dd.pi.degree.degree$space$plo,
                               dd.pi.degree.degree$space$phi,
                               dd.pi.degree.degree$time$plo,
                               dd.pi.degree.degree$time$phi)),
                   xs="r.degree",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=TRUE,
                   plot.x=FALSE,
                   treatments="all",
                   agg.col="GenusSpecies")
        mtext("Space", 3, line=1, cex=1.2)
        mtext(ylabs[1], 2,
              line=4, cex=1.2)
        ## time degree-degree
        plot.panel(new.dd=dd.pi.degree.degree$time,
                   dats=pol.specs[["time"]],
                   y1=ys[1],
                   y2= range(c(dd.pi.degree.degree$space$plo,
                               dd.pi.degree.degree$space$phi,
                               dd.pi.degree.degree$time$plo,
                               dd.pi.degree.degree$time$phi)),
                   xs="r.degree",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=FALSE,
                   plot.x=FALSE,
                   treatments="all",
                   agg.col="GenusSpecies")
        mtext("Time", 3, line=1, cex=1.2)

        ## space body size
        ## plot.panel(new.dd=dd.pi.body$space,
        ##            dats=pol.specs[["space"]],
        ##            y1=ys,
        ##            y2= range(c(dd.pi.body$space$plo,
        ##                        dd.pi.body$space$phi)),
        ##            xs="BodyLength",
        ##            col.fill=col.white,
        ##            col.lines=col.white,
        ##            col.points=col.fill,
        ##            plot.y=FALSE,
        ##            plot.x=FALSE,
        ##            treatments="all",
        ##            dec=0,
        ##            agg.col="GenusSpecies")

        ## space between-degree
        plot.panel(new.dd=dd.pi.degree.between$space,
                   dats=pol.specs[["space"]],
                   y1=ys[2],
                   y2= range(c(dd.pi.degree.between$space$plo,
                               dd.pi.degree.between$space$phi,
                               dd.pi.degree.between$time$plo,
                               dd.pi.degree.between$time$phi)),
                   xs="r.degree",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=TRUE,
                   plot.x=TRUE,
                   treatments="all",
                   agg.col="GenusSpecies")

        mtext(ylabs[2], 2,
              line=4, cex=1.2)
        mtext("Floral diet breadth", 1, line=3, cex=1.2)

        ## time between-degree
        plot.panel(new.dd=dd.pi.degree.between$time,
                   dats=pol.specs[["time"]],
                   y1=ys[2],
                   y2= range(c(dd.pi.degree.between$space$plo,
                               dd.pi.degree.between$space$phi,
                               dd.pi.degree.between$time$plo,
                               dd.pi.degree.between$time$phi)),
                   xs="r.degree",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=FALSE,
                   plot.x=TRUE,
                   treatments="all",
                   agg.col="GenusSpecies")

        mtext("Floral diet breadth", 1, line=3, cex=1.2)

        ## time bosy size
        ## plot.panel(new.dd=dd.pi.body.time,
        ##            dats=pol.specs[["time"]],
        ##            y1=ys,
        ##            y2= range(c(dd.pi.body$time$plo,
        ##                        dd.pi.body$time$phi)),
        ##            xs="BodyLength",
        ##            col.fill=col.white,
        ##            col.lines=col.white,
        ##            col.points=col.fill,
        ##            plot.y=FALSE,
        ##            plot.x=TRUE,
        ##            treatments="all",
        ##            dec=0,
        ##            agg.col="GenusSpecies")
    }
    path <- '../../../hedgerow_metacom_saved/occupancy/figures/ms'
    pdf.f(f, file=file.path(path,
                            sprintf("%s.pdf", "betweenness_degree")),
          width=5, height=4)
}





plot.panels.sites <- function(){
    f <- function(){
        col.lines <- "black"
        col.fill <- add.alpha(col.lines, alpha=0.2)
        col.white <- add.alpha("white", alpha=0)
        names(col.lines) <- names(col.fill) <- "all"
        ## layout(matrix(1:4, nrow=2, byrow=TRUE))
        layout(matrix(1, nrow=1, byrow=TRUE))
        par(oma=c(4.5, 7, 2.5, 1),
            mar=c(1, 0, 0, 1.5), cex.axis=1.2)

        ys <- c("k", "betweenness")
        ylabs <- c("Degree (z)", "Centrality")
        dd.pi.frd.degree <- makePredictData(frd.dd, ys[1],
                                            mods.site, time=FALSE)
        dd.pi.natArea.degree <- makePredictData(natArea.dd, ys[1],
                                                mods.site, time=FALSE)
        dd.pi.frd.between <- makePredictData(frd.dd, ys[2],
                                             mods.site, time=FALSE)
        dd.pi.natArea.between <- makePredictData(natArea.dd, ys[2],
                                                 mods.site, time=FALSE)
        ## space frd-degree
        plot.panel(new.dd=dd.pi.frd.degree$space,
                   dats=specs,
                   y1=ys[1],
                   y2= range(c(dd.pi.frd.degree$space$plo,
                               dd.pi.frd.degree$space$phi,
                               dd.pi.natArea.degree$space$plo,
                               dd.pi.natArea.degree$space$phi)),
                   xs="Div",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=TRUE,
                   plot.x=TRUE,
                   treatments="all",
                   agg.col="Site")
        mtext(ylabs[1], 2,
              line=4, cex=1.2)
        ## space nat area-degree
        ## plot.panel(new.dd=dd.pi.natArea.degree$space,
        ##            dats=specs,
        ##            y1=ys[1],
        ##            y2= range(c(dd.pi.frd.degree$space$plo,
        ##                        dd.pi.frd.degree$space$phi,
        ##                        dd.pi.natArea.degree$space$plo,
        ##                        dd.pi.natArea.degree$space$phi)),
        ##            xs="natArea",
        ##            col.fill=col.fill,
        ##            col.lines=col.lines,
        ##            plot.y=FALSE,
        ##            plot.x=FALSE,
        ##            treatments="all",
        ##            agg.col="Site")

        ## space frd-betweenness
        ## plot.panel(new.dd=dd.pi.frd.between$space,
        ##            dats=specs,
        ##            y1=ys[2],
        ##            y2= range(c(dd.pi.frd.between$space$plo,
        ##                        dd.pi.frd.between$space$phi,
        ##                        dd.pi.natArea.between$space$plo,
        ##                        dd.pi.natArea.between$space$phi)),
        ##            xs="Div",
        ##            col.fill=col.white,
        ##            col.lines=col.white,
        ##            col.points=col.fill,
        ##            plot.y=TRUE,
        ##            plot.x=TRUE,
        ##            treatments="all",
        ##            agg.col="Site")
        ## mtext(ylabs[2], 2,
        ##       line=4, cex=1.2)
        mtext("Floral diversity", 1, line=3, cex=1.2)
        ## space nat hab-betweenness
        ##     plot.panel(new.dd=dd.pi.natArea.between$space,
        ##                dats=specs,
        ##                y1=ys[2],
        ##                y2= range(c(dd.pi.natArea.between$space$plo,
        ##                            dd.pi.natArea.between$space$phi,
        ##                            dd.pi.frd.between$space$plo,
        ##                            dd.pi.frd.between$space$phi)),
        ##                xs="natArea",
        ##                col.fill=col.white,
        ##                col.lines=col.white,
        ##                col.points=col.fill,
        ##                plot.y=FALSE,
        ##                plot.x=TRUE,
        ##                treatments="all",
        ##                agg.col="Site")
        ##     mtext("Non-crop habitat \n proximity", 1, line=4, cex=1.2)
    }
    path <- '../../../hedgerow_metacom_saved/occupancy/figures/ms'
    pdf.f(f, file=file.path(path,
                            sprintf("%s.pdf", "betweenness_degree_sites")),
          width=4, height=3.5)
}
