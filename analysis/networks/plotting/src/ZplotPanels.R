
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

plot.panels.all <- function(){
    f <- function(){
        col.lines <- "black"
        col.fill <- add.alpha(col.lines, alpha=0.2)
        col.white <- add.alpha("white", alpha=0)
        names(col.lines) <- names(col.fill) <- "all"
        layout(matrix(1:(5*2), nrow=2, byrow=TRUE))
        par(oma=c(4, 7, 3, 1),
            mar=c(1, 0, 0, 1.5), cex.axis=1.2)

        ys <- c("degree", "betweenness")
        ylabs <- c("Degree", "Betweenness")

        dd.pi.degree.degree <- makePredictData(r.degree.dd, ys[1], pol.mods)
        dd.pi.degree.between <- makePredictData(r.degree.dd, ys[2],
                                                pol.mods)
        ## dd.pi.degree.close <- makePredictData(r.degree.dd, ys[3],
        ##                                       pol.mods)

        dd.pi.body.degree <- makePredictData(body.dd, ys[1], pol.mods)
        dd.pi.body.between <- makePredictData(body.dd, ys[2],
                                              pol.mods)
        ## dd.pi.body.close <- makePredictData(body.dd, ys[3],
        ##                                     pol.mods)

        dd.pi.frd.degree <- makePredictData(frd.dd, ys[1],
                                            mods.site, time=FALSE)
        dd.pi.natArea.degree <- makePredictData(natArea.dd, ys[1],
                                                mods.site, time=FALSE)


        dd.pi.frd.between <- makePredictData(frd.dd, ys[2],
                                             mods.site, time=FALSE)
        dd.pi.natArea.between <- makePredictData(natArea.dd, ys[2],
                                                 mods.site,
                                                 time=FALSE)

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
        ## mtext("Spatial: pollinators", 3, line=1, cex=1.2)
        mtext(ylabs[1], 2,
              line=4, cex=1.2)

        ## space body-degree
        plot.panel(new.dd=dd.pi.body.degree$space,
                   dats=pol.specs[["space"]],
                   y1=ys[1],
                   y2= range(c(dd.pi.body.degree$space$plo,
                               dd.pi.body.degree$space$phi,
                               dd.pi.body.degree$time$plo,
                               dd.pi.body.degree$time$phi)),
                   xs="MeanITD",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=FALSE,
                   plot.x=FALSE,
                   treatments="all",
                   agg.col="GenusSpecies")

        mtext("Spatial: pollinators", 3, line=1, cex=1.2, adj=-4)

        ## time degree-degree
        plot.panel(new.dd=dd.pi.degree.degree$time,
                   dats=pol.specs[["time"]],
                   y1=ys[1],
                   y2= range(c(dd.pi.degree.degree$space$plo,
                               dd.pi.degree.degree$space$phi,
                               dd.pi.degree.degree$time$plo,
                               dd.pi.degree.degree$time$phi, 0)),
                   xs="r.degree",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=FALSE,
                   plot.x=FALSE,
                   treatments="all",
                   agg.col="GenusSpecies")
        ## mtext("Temporal: pollinators", 3, line=1, cex=1.2)


        ## time degree-degree
        plot.panel(new.dd=dd.pi.body.degree$time,
                   dats=pol.specs[["time"]],
                   y1=ys[1],
                   y2= range(c(dd.pi.body.degree$space$plo,
                               dd.pi.body.degree$space$phi,
                               dd.pi.body.degree$time$plo,
                               dd.pi.body.degree$time$phi, 0)),
                   xs="MeanITD",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=FALSE,
                   plot.x=FALSE,
                   treatments="all",
                   agg.col="GenusSpecies")
        mtext("Temporal: pollinators", 3, line=1, cex=1.2, adj=-1200)



        ## space frd-degree
        plot.panel(new.dd=dd.pi.frd.degree$space,
                   dats=specs,
                   y1=ys[1],
                   y2= range(c(dd.pi.frd.degree$space$plo,
                               dd.pi.frd.degree$space$phi,
                               dd.pi.natArea.degree$space$plo,
                               dd.pi.natArea.degree$space$phi, 0)),
                   xs="Div",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=FALSE,
                   plot.x=FALSE,
                   treatments="all",
                   agg.col="Site")
        ## mtext(ylabs[1], 2,
        ##       line=4, cex=1.2)
        mtext("Spatial: sites", 3, line=1, cex=1.2)

################################################
        ##next row betweenness
        ## space close-degree
        plot.panel(new.dd=dd.pi.degree.between$space,
                   dats=pol.specs[["space"]],
                   y1=ys[2],
                   y2= range(c(dd.pi.degree.between$space$plo,
                               dd.pi.degree.between$space$phi,
                               dd.pi.degree.between$time$plo,
                               dd.pi.degree.between$time$phi, 0)),
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

        ## space between-body
        plot.panel(new.dd=dd.pi.body.between$space,
                   dats=pol.specs[["space"]],
                   y1=ys[2],
                   y2= range(c(dd.pi.body.between$space$plo,
                               dd.pi.body.between$space$phi,
                               dd.pi.body.between$time$plo,
                               dd.pi.body.between$time$phi, 0)),
                   xs="MeanITD",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   col.points=col.fill,
                   plot.y=FALSE,
                   plot.x=TRUE,
                   treatments="all",
                   agg.col="GenusSpecies")

               mtext("Body size", 1, line=3, cex=1.2)

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

        ## time between-degree
        plot.panel(new.dd=dd.pi.body.between$time,
                   dats=pol.specs[["time"]],
                   y1=ys[2],
                   y2= range(c(dd.pi.body.between$space$plo,
                               dd.pi.body.between$space$phi,
                               dd.pi.body.between$time$plo,
                               dd.pi.body.between$time$phi)),
                   xs="MeanITD",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=FALSE,
                   plot.x=TRUE,
                   treatments="all",
                   agg.col="GenusSpecies")

                       mtext("Body size", 1, line=3, cex=1.2)


        plot.panel(new.dd=dd.pi.frd.between$space,
                   dats=specs,
                   y1=ys[2],
                   y2= range(c(dd.pi.frd.between$space$plo,
                               dd.pi.frd.between$space$phi,
                               dd.pi.natArea.between$space$plo,
                               dd.pi.natArea.between$space$phi, 0)),
                   xs="Div",
                   col.fill=col.white,
                   col.lines=col.white,
                   col.points=col.fill,
                   plot.y=FALSE,
                   plot.x=TRUE,
                   treatments="all",
                   agg.col="Site")

                mtext("Floral diversity", 1, line=3, cex=1.2)
        ## mtext(ylabs[2], 2,
        ##       line=4, cex=1.2)

## ################################################
##         ##next row
##         ## space close-degree
##         plot.panel(new.dd=dd.pi.degree.close$space,
##                    dats=pol.specs[["space"]],
##                    y1=ys[3],
##                    y2= range(c(dd.pi.degree.close$space$plo,
##                                dd.pi.degree.close$space$phi,
##                                dd.pi.degree.close$time$plo,
##                                dd.pi.degree.close$time$phi)),
##                    xs="r.degree",
##                    col.fill=col.fill,
##                    col.lines=col.lines,
##                    plot.y=TRUE,
##                    plot.x=TRUE,
##                    treatments="all",
##                    agg.col="GenusSpecies")

##         mtext(ylabs[3], 2,
##               line=4, cex=1.2)
##         mtext("Floral diet breadth", 1, line=3, cex=1.2)

##         ## space close-body
##         plot.panel(new.dd=dd.pi.body.close$space,
##                    dats=pol.specs[["space"]],
##                    y1=ys[3],
##                    y2= range(c(dd.pi.body.close$space$plo,
##                                dd.pi.body.close$space$phi,
##                                dd.pi.body.close$time$plo,
##                                dd.pi.body.close$time$phi)),
##                    xs="MeanITD",
##                    col.fill=col.fill,
##                    col.lines=col.lines,
##                    col.points=col.fill,
##                    plot.y=FALSE,
##                    plot.x=TRUE,
##                    treatments="all",
##                    agg.col="GenusSpecies")

##         mtext("Body size", 1, line=3, cex=1.2)

##         ## time close-degree
##         plot.panel(new.dd=dd.pi.degree.close$time,
##                    dats=pol.specs[["time"]],
##                    y1=ys[3],
##                    y2= range(c(dd.pi.degree.close$space$plo,
##                                dd.pi.degree.close$space$phi,
##                                dd.pi.degree.close$time$plo,
##                                dd.pi.degree.close$time$phi)),
##                    xs="r.degree",
##                    col.fill=col.fill,
##                    col.lines=col.lines,
##                    plot.y=FALSE,
##                    plot.x=TRUE,
##                    treatments="all",
##                    agg.col="GenusSpecies")

##         mtext("Floral diet breadth", 1, line=3, cex=1.2)



##         ## time close-degree
##         plot.panel(new.dd=dd.pi.body.close$time,
##                    dats=pol.specs[["time"]],
##                    y1=ys[3],
##                    y2= range(c(dd.pi.body.close$space$plo,
##                                dd.pi.body.close$space$phi,
##                                dd.pi.body.close$time$plo,
##                                dd.pi.body.close$time$phi)),
##                    xs="MeanITD",
##                    col.fill=col.fill,
##                    col.lines=col.lines,
##                    plot.y=FALSE,
##                    plot.x=TRUE,
##                    treatments="all",
##                    agg.col="GenusSpecies")

##         mtext("Body size", 1, line=3, cex=1.2)


##         plot.panel(new.dd=dd.pi.frd.close$space,
##                    dats=specs,
##                    y1=ys[3],
##                    y2= range(c(dd.pi.frd.close$space$plo,
##                                dd.pi.frd.close$space$phi,
##                                dd.pi.natArea.close$space$plo,
##                                dd.pi.natArea.close$space$phi)),
##                    xs="Div",
##                    col.fill=col.white,
##                    col.lines=col.white,
##                    col.points=col.fill,
##                    plot.y=FALSE,
##                    plot.x=TRUE,
##                    treatments="all",
##                    agg.col="Site")
##         ## mtext(ylabs[2], 2,
##         ##       line=4, cex=1.2)
##         mtext("Floral diversity", 1, line=3, cex=1.2)

    }
    path <- '../../../hedgerow_metacom_saved/networks/figures'
    pdf.f(f, file=file.path(path,
                            sprintf("all_sig_drop_li_ht%s.pdf", drop.li.ht)),
          width=11, height=5)
}




