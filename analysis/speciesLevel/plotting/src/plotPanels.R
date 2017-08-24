plot.panels <- function(){
    f <- function(){
        col.lines <- "black"
        col.fill <- add.alpha(col.lines, alpha=0.2)
        col.white <- add.alpha("white", alpha=0)
        names(col.lines) <- names(col.fill) <- "all"
        layout(matrix(1:4, nrow=2, byrow=TRUE))
        par(oma=c(4, 6.5, 0.5, 1),
            mar=c(1, 0, 0, 1.5), cex.axis=1.2)

        ys <- c("k")
        xvar <- c("r.degree", "BodyLength")
        ylabs <- c("K")

        dd.degree <- cbind(r.degree.dd, 0)
        colnames(dd.degree) <- c(colnames(r.degree.dd), ys)
        dd.pi.degree.space <- predict.int(mod= mods[["space"]][[ys]],
                                          dd=dd.degree,
                                          y=ys,
                                          family="gaussian")
        dd.pi.degree.time <- predict.int(mod= mods[["time"]][[ys]],
                                         dd=dd.degree,
                                         y=ys,
                                         family="gaussian")

        dd.body <- cbind(body.dd, 0)
        colnames(dd.body) <- c(colnames(body.dd), ys)
        dd.pi.body.space <- predict.int(mod= mods[["space"]][[ys]],
                                        dd=dd.body,
                                        y=ys,
                                        family="gaussian")
        dd.pi.body.time <- predict.int(mod= mods[["time"]][[ys]],
                                       dd=dd.body,
                                       y=ys,
                                       family="gaussian")


        ## space degree
        plot.panel(new.dd=dd.pi.degree.space,
                   dats=all.specs[["space"]],
                   y1=ys,
                   y2= range(c(dd.pi.degree.space$plo,
                               dd.pi.degree.space$phi)),
                   xs="r.degree",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=TRUE,
                   plot.x=FALSE,
                   treatments="all",
                   agg.col="GenusSpecies")

        mtext("Space", 2, line=3.5, cex=1.2)
        mtext("K", 2,
              line=2, cex=1.2)

        ## space body size
        plot.panel(new.dd=dd.pi.body.space,
                   dats=all.specs[["space"]],
                   y1=ys,
                   y2= range(c(dd.pi.body.space$plo,
                               dd.pi.body.space$phi)),
                   xs="BodyLength",
                   col.fill=col.white,
                   col.lines=col.white,
                   col.points=col.fill,
                   plot.y=FALSE,
                   plot.x=FALSE,
                   treatments="all",
                   dec=0,
                   agg.col="GenusSpecies")

        ## time degree
        plot.panel(new.dd=dd.pi.degree.time,
                   dats=all.specs[["time"]],
                   y1=ys,
                   y2= range(c(dd.pi.degree.time$plo,
                               dd.pi.degree.time$phi)),
                   xs="r.degree",
                   col.fill=col.fill,
                   col.lines=col.lines,
                   plot.y=TRUE,
                   plot.x=TRUE,
                   treatments="all",
                   agg.col="GenusSpecies")

        mtext("Time", 2, line=3.5, cex=1.2)

        mtext("Floral Degree", 1, line=3, cex=1.2)
        mtext("K", 2,
              line=2, cex=1.2)


        ## degree plants
        plot.panel(new.dd=dd.pi.body.time,
                   dats=all.specs[["time"]],
                   y1=ys,
                   y2= range(c(dd.pi.body.time$plo,
                               dd.pi.body.time$phi)),
                   xs="BodyLength",
                   col.fill=col.white,
                   col.lines=col.white,
                   col.points=col.fill,
                   plot.y=FALSE,
                   plot.x=TRUE,
                   treatments="all",
                   dec=0,
                   agg.col="GenusSpecies")
        mtext("Body size", 1, line=3, cex=1.2)
    }
    path <- 'figures'
    pdf.f(f, file=file.path(path,
                            sprintf("%s.pdf", "degree_bodysize")),
          width=5, height=4)
}

