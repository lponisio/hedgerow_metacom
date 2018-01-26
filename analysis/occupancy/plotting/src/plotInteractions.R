
plotHRInteractions <- function(){
    layout(matrix(c(1,1,1, 2:7), nrow=3, byrow=TRUE),
           heights=c(0.25,1,1,1,1,1,1))

    par(oma=c(3, 7, 2.5, 1),
        mar=c(5, 0, 0, 1.5), cex.axis=1.5)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", bty="n", ylab="", xlab="")
    legend("center", legend=c("min", "2.5%", "25%", "median", "75%", "97.5%", "max"),
           pch=16, col=na.omit(cols[match(probs, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1))]),
           bty="n", cex=1.5, ncol=7)

    ## 1
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.fra'] * quantiles$fra[i] +
                        means['mu.phi.hr.area.fra'] * x * quantiles$fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Floral diversity (a)", bty="n", cex=1.2)

    ## 2
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="", yaxt="n")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.k'] * quantiles$k[i] +
                        means['mu.phi.hr.area.k'] * x * quantiles$k[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Diet breadth (b)", bty="n", cex=1.2)
    ## 3
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         yaxt="n", ylab="", xlab="")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.B'] * quantiles$B[i] +
                        means['mu.phi.hr.area.B'] * x * quantiles$B[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Body size (c)", bty="n", cex=1.2)

    ## 4
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Colonization", 2, line=4, cex=1.3)
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.fra'] * quantiles$fra[i] +
                        means['mu.gam.hr.area.fra'] * x * quantiles$fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topleft", legend="Floral diversity (a)", bty="n", cex=1.2)

    ## 5
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="", yaxt="n")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.k'] * quantiles$k[i] +
                        means['mu.gam.hr.area.k'] * x * quantiles$k[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topleft", legend="Diet breadth (b)", bty="n", cex=1.2)

    ## 6
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         yaxt="n", ylab="", xlab="")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.B'] * quantiles$B[i] +
                        means['mu.gam.hr.area.B'] * x * quantiles$B[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topleft", legend="Body size (c)", bty="n", cex=1.2)
}




plotHRPersistence <- function(){
    layout(matrix(c(1,1,1, 2:4), nrow=2, byrow=TRUE), heights=c(0.25,1,1,1))
    par(oma=c(6, 7, 2.5, 1),
        mar=c(1, 0, 0, 1.5), cex.axis=1.5)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", bty="n", ylab="", xlab="")
    legend("center", legend=c("min", "5%", "25%", "median", "75%", "95%", "max"),
           pch=16, col=na.omit(cols[match(probs, c(0, 0.025, 0.25, 0.5, 0.75, 0.95, 1))]),
           bty="n", cex=1.5, ncol=7)

    ## 1
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.fra'] * quantiles$fra[i] +
                        means['mu.phi.hr.area.fra'] * x * quantiles$fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="Floral diversity (a)", bty="n", cex=1.2)

    ## 2
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="", yaxt="n")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.k'] * quantiles$k[i] +
                        means['mu.phi.hr.area.k'] * x * quantiles$k[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Diet breadth (b)", bty="n", cex=1.2)
    ## 3
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         yaxt="n", ylab="", xlab="")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.B'] * quantiles$B[i] +
                        means['mu.phi.hr.area.B'] * x * quantiles$B[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Body size (c)", bty="n", cex=1.2)
}





plotHRPerstTurnOccCol <- function(){
    layout(matrix(c(1,1, 2:5), nrow=3, byrow=TRUE),
           heights=c(0.75,1,1,1,1))

    par(oma=c(2, 7, 2, 1),
        mar=c(4.5, 5, 1, 1.5), cex.axis=1.5)

    kept.traits <- all.traits$r.degree[all.traits$GenusSpecies %in%
                                       names(model.input$data$k)]
    nbreaks <- 10
    h <- hist(kept.traits, breaks=nbreaks, plot=FALSE)
    cols <- rev(viridis(length(h$density)))
    plot(h, col=cols,
         xlab="", main="", ylab="", las=1)
    abline(v=mean(kept.traits), lty=2)
    ## legend("topright", legend=c("min", "2.5%", "25%", "median", "75%", "97.5%", "max"),
    ##        pch=16, col=cols[match(c(0.000, 0.025, 0.250, 0.500, 0.750, 0.9, 1.000), probs)],
    ##        bty="n", cex=1.5, ncol=7)

    quantiles$k <- (h$mids - mean(kept.traits))/sd(kept.traits)

    mtext("Frequency", 2, line=4, cex=1.3)
    mtext("Floral diet breadth", 1, line=3, cex=1.3)
    legend("topright", legend="(a)", bty="n", cex=1.2)

    ## 1 persistence
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    legend("bottomleft", legend="(b)", bty="n", cex=1.2)
    ## mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(quantiles$k)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.k'] * quantiles$k[i] +
                        means['mu.phi.hr.area.k'] * x * quantiles$k[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 2 colonization
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    legend("topleft", legend="(c)", bty="n", cex=1.2)
    mtext("Colonization", 2, line=4, cex=1.3)

    for(i in 1:length(quantiles$k)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.k'] * quantiles$k[i] +
                        means['mu.gam.hr.area.k'] * x * quantiles$k[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 3 occupancy
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         ylab="", xlab="", las=1)
    legend("topleft", legend="(d)", bty="n", cex=1.2)
    mtext("Occupancy", 2, line=4, cex=1.3)
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(quantiles$k)){
        curve((inv.logit(means['mu.gam.0'] +
                         means['mu.gam.hr.area'] * x +
                         means['mu.gam.k'] * quantiles$k[i] +
                         means['mu.gam.hr.area.k'] * x *
                         quantiles$k[i]))/
              (1 + inv.logit(means['mu.gam.0'] +
                             means['mu.gam.hr.area'] * x +
                             means['mu.gam.k'] * quantiles$k[i] +
                             means['mu.gam.hr.area.k'] * x *
                             quantiles$k[i]) -
               inv.logit(means['mu.phi.0'] +
                         means['mu.phi.hr.area'] * x +
                         means['mu.phi.k'] * quantiles$k[i] +
                         means['mu.phi.hr.area.k'] * x * quantiles$k[i])),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 4 turnover
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea),
         las=1, ylab="", xlab="")
    legend("bottomleft", legend="(e)", bty="n", cex=1.2)
    mtext("Turnover", 2, line=4, cex=1.3)
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(quantiles$k)){
        curve(2*inv.logit(means['mu.gam.0'] +
                          means['mu.gam.hr.area'] * x +
                          means['mu.gam.k'] * quantiles$k[i] +
                          means['mu.gam.hr.area.k'] * x * quantiles$k[i])*
              (1 - inv.logit(means['mu.phi.0'] +
                             means['mu.phi.hr.area'] * x +
                             means['mu.phi.k'] * quantiles$k[i] +
                             means['mu.phi.hr.area.k'] * x *
                             quantiles$k[i]))/
              (1+inv.logit(means['mu.gam.0'] +
                           means['mu.gam.hr.area'] * x +
                           means['mu.gam.k'] * quantiles$k[i] +
                           means['mu.gam.hr.area.k'] * x *
                           quantiles$k[i]) -
               inv.logit(means['mu.phi.0'] +
                         means['mu.phi.hr.area'] * x +
                         means['mu.phi.k'] * quantiles$k[i] +
                         means['mu.phi.hr.area.k'] * x * quantiles$k[i])) ,
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
}



