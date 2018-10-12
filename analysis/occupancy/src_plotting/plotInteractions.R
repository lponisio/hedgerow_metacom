f.plotInteractionsHRRemnant.k <- function(){
    plotInteractionsHRRemnant(means=means,
                              all.traits=all.traits,
                              all.traits.col.name='r.degree',
                              model.input=model.input,
                              model.input.trait.col.name='k',
                              param.name='k',
                              nbreaks=30,
                              hist.xlab='Floral diet breadth',
                              legend.loc=c(rep("bottomleft", 2),
                                           rep("topright", 2),
                                           rep("bottomright", 2)))
}


plotInteractionsHRRemnant <- function(means,
                                      all.traits,
                                      all.traits.col.name,
                                      model.input,
                                      model.input.trait.col.name,
                                      param.name,
                                      nbreaks,
                                      hist.xlab,
                                      legend.loc){
    layout(matrix(c(1,1, 2:7), nrow=4, byrow=TRUE),
           heights=c(0.75,1,1,1))
    par(oma=c(5, 5, 2, 1),
        mar=c(4, 2, 1, 1.5), cex.axis=1.5)
    kept.traits <- all.traits[, all.traits.col.name][all.traits$GenusSpecies %in%
                                                     names(model.input$data[[model.input.trait.col.name]])]
    h <- hist(kept.traits, breaks=nbreaks, plot=FALSE)
    cols <- rev(viridis(length(h$density)))
    plot(h, col=cols,
         xlab="", main="", ylab="", las=1)
    abline(v=mean(kept.traits), lty=2, col="red", lwd=3)
    quantiles.trait <- (h$mids - mean(kept.traits))/sd(kept.traits)
    mtext("Frequency", 2, line=4, cex=1.3)
    mtext(hist.xlab, 1, line=3, cex=1)

    ## 1 persistence
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    legend(legend.loc[1], legend="(a)", bty="n", cex=1)
    for(i in 1:length(quantiles.trait)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means[paste0('phi.', param.name)] * quantiles.trait[i] +
                        means[paste0('phi.hr.area.', param.name)] * x * quantiles.trait[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 2 persistence
    ## interactions of floral resources and nat habitat proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural), las=1,
         ylab="", xlab="", yaxt="n")
    legend(legend.loc[2], legend="(b)", bty="n", cex=1)
    for(i in 1:length(quantiles.trait)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.nat.area'] * x +
                        means[paste0('phi.', param.name)] * quantiles.trait[i] +
                        means[paste0('phi.nat.area.', param.name)] * x * quantiles.trait[i]),
              from=range(model.input$data$natural)[1],
              to=range(model.input$data$natural)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 3 colonization
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    legend(legend.loc[3], legend="(c)", bty="n", cex=1)
    mtext("Colonization", 2, line=4, cex=1.3)
    ## mtext("Hedgerow proximity-weighted \n area", 1, line=5, cex=1.3)
    for(i in 1:length(quantiles.trait)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means[paste0('gam.', param.name)] * quantiles.trait[i] +
                        means[paste0('gam.hr.area.', param.name)] * x * quantiles.trait[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 4 colonization
    ## interactions of floral resources and natural habitat
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural), las=1,
         ylab="", xlab="", yaxt="n")
    legend(legend.loc[4], legend="(d)", bty="n", cex=1)

    ## mtext("Remnant proximity-weighted \n area", 1, line=5, cex=1.3)
    for(i in 1:length(quantiles.trait)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.nat.area'] * x +
                        means[paste0('gam.', param.name)] * quantiles.trait[i] +
                        means[paste0('gam.nat.area.', param.name)] * x * quantiles.trait[i]),
              from=range(model.input$data$natural)[1],
              to=range(model.input$data$natural)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }


    ## 5 occupancy hedgerow*k
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         ylab="", xlab="", las=1)
    legend(legend.loc[5], legend="(e)", bty="n", cex=1)
    mtext("Occupancy", 2, line=4, cex=1.3)

    for(i in 1:length(quantiles.trait)){
        curve((inv.logit(means['mu.gam.0'] +
                         means['mu.gam.hr.area'] * x +
                         means[paste0('gam.', param.name)] * quantiles.trait[i] +
                         means[paste0('gam.hr.area.', param.name)] * x *
                         quantiles.trait[i]))/
              (1 + inv.logit(means['mu.gam.0'] +
                             means['mu.gam.hr.area'] * x +
                             means[paste0('gam.', param.name)] * quantiles.trait[i] +
                             means[paste0('gam.hr.area.', param.name)] * x *
                             quantiles.trait[i]) -
               inv.logit(means['mu.phi.0'] +
                         means['mu.phi.hr.area'] * x +
                         means[paste0('phi.', param.name)] * quantiles.trait[i] +
                         means[paste0('phi.hr.area.', param.name)] * x * quantiles.trait[i])),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }


    mtext("Hedgerow \n proximity-weighted area", 1, line=5.5, cex=1)

    ## 6 occupancy nat hab*k
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$natural),
         ylab="", xlab="", las=1, yaxt="n")
    legend(legend.loc[6], legend="(f)", bty="n", cex=1)
    ## mtext("Occupancy", 2, line=4, cex=1.3)

    for(i in 1:length(quantiles.trait)){
        curve((inv.logit(means['mu.gam.0'] +
                         means['mu.gam.nat.area'] * x +
                         means[paste0('gam.', param.name)] * quantiles.trait[i] +
                         means[paste0('gam.nat.area.', param.name)] * x *
                         quantiles.trait[i]))/
              (1 + inv.logit(means['mu.gam.0'] +
                             means['mu.gam.nat.area'] * x +
                             means[paste0('gam.', param.name)] * quantiles.trait[i] +
                             means[paste0('gam.nat.area.', param.name)] * x *
                             quantiles.trait[i]) -
               inv.logit(means['mu.phi.0'] +
                         means['mu.phi.nat.area'] * x +
                         means[paste0('phi.', param.name)] * quantiles.trait[i] +
                         means[paste0('phi.nat.area.', param.name)] * x * quantiles.trait[i])),
              from=range(model.input$data$natural)[1],
              to=range(model.input$data$natural)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    mtext("Remnant habitat \n proximity-weighted area", 1, line=5.5, cex=1)
}

## *********************************************************************************

f.plotInteractionsFloralDiv <- function(){
    plotInteractionsFloralDiv(plot.turnover=FALSE,
                              plot.remnant=FALSE)
}

plotInteractionsFloralDiv <- function(plot.turnover,
                                      plot.remnant){
    if(plot.turnover){
        layout(matrix(c(1,1, 2:5), nrow=3, byrow=TRUE),
               heights=c(0.75,1,1,1,1))
    } else if(plot.remnant){
        layout(matrix(c(1,1, 2:7), nrow=4, byrow=TRUE),
               heights=c(0.75,1,1,1))
    } else{
        layout(matrix(c(1,1, 2,2, 3,3, 4,4), nrow=4, byrow=TRUE),
               heights=c(0.75,1,1,1))
    }
    par(oma=c(5, 5, 2, 1),
        mar=c(4, 2, 1, 1.5), cex.axis=1.5)
    nbreaks <- 30

    h <- hist(model.input$data$fra, breaks=nbreaks, plot=FALSE)
    cols <- rev(viridis(length(h$density)))
    plot(h, col=cols,
         xlab="", main="", ylab="", las=1)
    abline(v=mean(model.input$data$fra), lty=2, col="red", lwd=3)
    quantiles.fra <- (h$mids - mean(model.input$data$fra))/sd(model.input$data$fra)
    mtext("Frequency", 2, line=4, cex=1.3)
    mtext("Floral diversity", 1, line=3, cex=1)

    ## 1 persistence
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    legend("topleft", legend="(a)", bty="n", cex=1.2)
    for(i in 1:length(quantiles.fra)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.fra'] * quantiles.fra[i] +
                        means['phi.hr.area.fra'] * x * quantiles.fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    if(plot.remnant){
        ## 2 persistence
        ## interactions of floral resources and nat habitat proximity
        plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural), las=1,
             ylab="", xlab="", yaxt="n")
        legend("topright", legend="(b)", bty="n", cex=1.2)
        for(i in 1:length(quantiles.fra)){
            curve(inv.logit(means['mu.phi.0'] +
                            means['mu.phi.nat.area'] * x +
                            means['mu.phi.fra'] * quantiles.fra[i] +
                            means['phi.nat.area.fra'] * x * quantiles.fra[i]),
                  from=range(model.input$data$natural)[1],
                  to=range(model.input$data$natural)[2],
                  col=cols[i],
                  lwd=2,
                  add=TRUE)
        }
    }

    ## 3 colonization
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    legend("topleft", legend="(b)", bty="n", cex=1.2)
    mtext("Colonization", 2, line=4, cex=1.3)
    ## mtext("Hedgerow proximity-weighted \n area", 1, line=5, cex=1.3)
    for(i in 1:length(quantiles.fra)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.fra'] * quantiles.fra[i] +
                        means['gam.hr.area.fra'] * x * quantiles.fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    if(plot.remnant){
        ## 4 colonization
        ## interactions of floral resources and natural habitat
        plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural), las=1,
             ylab="", xlab="", yaxt="n")
        legend("topleft", legend="(d)", bty="n", cex=1.2)

        ## mtext("Remnant proximity-weighted \n area", 1, line=5, cex=1.3)
        for(i in 1:length(quantiles.fra)){
            curve(inv.logit(means['mu.gam.0'] +
                            means['mu.gam.nat.area'] * x +
                            means['mu.gam.fra'] * quantiles.fra[i] +
                            means['gam.nat.area.fra'] * x * quantiles.fra[i]),
                  from=range(model.input$data$natural)[1],
                  to=range(model.input$data$natural)[2],
                  col=cols[i],
                  lwd=2,
                  add=TRUE)
        }
    }

    ## 5 occupancy hedgerow*fra
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         ylab="", xlab="", las=1)
    legend("topleft", legend="(c)", bty="n", cex=1.2)
    mtext("Occupancy", 2, line=4, cex=1.3)

    for(i in 1:length(quantiles.fra)){
        curve((inv.logit(means['mu.gam.0'] +
                         means['mu.gam.hr.area'] * x +
                         means['mu.gam.fra'] * quantiles.fra[i] +
                         means['gam.hr.area.fra'] * x *
                         quantiles.fra[i]))/
              (1 + inv.logit(means['mu.gam.0'] +
                             means['mu.gam.hr.area'] * x +
                             means['mu.gam.fra'] * quantiles.fra[i] +
                             means['gam.hr.area.fra'] * x *
                             quantiles.fra[i]) -
               inv.logit(means['mu.phi.0'] +
                         means['mu.phi.hr.area'] * x +
                         means['mu.phi.fra'] * quantiles.fra[i] +
                         means['phi.hr.area.fra'] * x * quantiles.fra[i])),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
        mtext("Hedgerow \n proximity-weighted area", 1, line=5.5, cex=1.2)
    }
    if(plot.remnant){
        ## 6 occupancy nat hab*fra
        plot(NA, ylim=c(0, 1),
             xlim=range(model.input$data$natural),
             ylab="", xlab="", las=1, yaxt="n")
        legend("topleft", legend="(f)", bty="n", cex=1.2)
        ## mtext("Occupancy", 2, line=4, cex=1.3)

        for(i in 1:length(quantiles.fra)){
            curve((inv.logit(means['mu.gam.0'] +
                             means['mu.gam.nat.area'] * x +
                             means['mu.gam.fra'] * quantiles.fra[i] +
                             means['gam.nat.area.fra'] * x *
                             quantiles.fra[i]))/
                  (1 + inv.logit(means['mu.gam.0'] +
                                 means['mu.gam.nat.area'] * x +
                                 means['mu.gam.fra'] * quantiles.fra[i] +
                                 means['gam.nat.area.fra'] * x *
                                 quantiles.fra[i]) -
                   inv.logit(means['mu.phi.0'] +
                             means['mu.phi.nat.area'] * x +
                             means['mu.phi.fra'] * quantiles.fra[i] +
                             means['phi.nat.area.fra'] * x * quantiles.fra[i])),
                  from=range(model.input$data$natural)[1],
                  to=range(model.input$data$natural)[2],
                  col=cols[i],
                  lwd=2,
                  add=TRUE)
        }
    }

    if(plot.turnover){
        ## 7 turnover hedgerow*fra
        plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea),
             las=1, ylab="", xlab="")
        legend("topleft", legend="(g)", bty="n", cex=1.2)
        mtext("Turnover", 2, line=4, cex=1.3)
        mtext("Hedgerow \n proximity-weighted area", 1, line=5.5, cex=1.2)

        for(i in 1:length(quantiles.fra)){
            curve(2*inv.logit(means['mu.gam.0'] +
                              means['mu.gam.hr.area'] * x +
                              means['mu.gam.fra'] * quantiles.fra[i] +
                              means['gam.hr.area.fra'] * x * quantiles.fra[i])*
                  (1 - inv.logit(means['mu.phi.0'] +
                                 means['mu.phi.hr.area'] * x +
                                 means['mu.phi.fra'] * quantiles.fra[i] +
                                 means['phi.hr.area.fra'] * x *
                                 quantiles.fra[i]))/
                  (1+inv.logit(means['mu.gam.0'] +
                               means['mu.gam.hr.area'] * x +
                               means['mu.gam.fra'] * quantiles.fra[i] +
                               means['gam.hr.area.fra'] * x *
                               quantiles.fra[i]) -
                   inv.logit(means['mu.phi.0'] +
                             means['mu.phi.hr.area'] * x +
                             means['mu.phi.fra'] * quantiles.fra[i] +
                             means['phi.hr.area.fra'] * x * quantiles.fra[i])) ,
                  from=range(model.input$data$HRarea)[1],
                  to=range(model.input$data$HRarea)[2],
                  col=cols[i],
                  lwd=2,
                  add=TRUE)
            mtext("Hedgerow \n proximity-weighted area", 1, line=5.5, cex=1.2)
        }
    }
    if(plot.remnant){

        ## 8 turnover nat*fra
        plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural),
             las=1, ylab="", xlab="", yaxt="n")
        legend("topleft", legend="(h)", bty="n", cex=1.2)
        ## mtext("Turnover", 2, line=4, cex=1.3)
        mtext("Remnant \n proximity-weighted area", 1, line=5.5, cex=1.2)

        for(i in 1:length(quantiles.fra)){
            curve(2*inv.logit(means['mu.gam.0'] +
                              means['mu.gam.nat.area'] * x +
                              means['mu.gam.fra'] * quantiles.fra[i] +
                              means['gam.nat.area.fra'] * x * quantiles.fra[i])*
                  (1 - inv.logit(means['mu.phi.0'] +
                                 means['mu.phi.nat.area'] * x +
                                 means['mu.phi.fra'] * quantiles.fra[i] +
                                 means['phi.nat.area.fra'] * x *
                                 quantiles.fra[i]))/
                  (1+inv.logit(means['mu.gam.0'] +
                               means['mu.gam.nat.area'] * x +
                               means['mu.gam.fra'] * quantiles.fra[i] +
                               means['gam.nat.area.fra'] * x *
                               quantiles.fra[i]) -
                   inv.logit(means['mu.phi.0'] +
                             means['mu.phi.nat.area'] * x +
                             means['mu.phi.fra'] * quantiles.fra[i] +
                             means['phi.nat.area.fra'] * x * quantiles.fra[i])) ,
                  from=range(model.input$data$natural)[1],
                  to=range(model.input$data$natural)[2],
                  col=cols[i],
                  lwd=2,
                  add=TRUE)
            mtext("Remnant \n proximity-weighted area", 1, line=5.5, cex=1.2)
        }

    }
}


## *********************************************************************************


plotInteractionsB <- function(){

    layout(matrix(c(1,1, 2,2, 3,3, 4,4), nrow=4, byrow=TRUE),
           heights=c(0.75,1,1,1))
    par(oma=c(5, 5, 2, 1),
        mar=c(4, 2, 1, 1.5), cex.axis=1.5)
    kept.traits <- all.traits$MeanITD[all.traits$GenusSpecies %in%
                                      names(model.input$data$B)]
    nbreaks <- 30
    h <- hist(kept.traits, breaks=nbreaks, plot=FALSE)
    cols <- rev(viridis(length(h$density)))
    plot(h, col=cols,
         xlab="", main="", ylab="", las=1)
    abline(v=mean(kept.traits), lty=2, col="red", lwd=3)

    quantiles.B <- (h$mids - mean(kept.traits))/sd(kept.traits)

    mtext("Frequency", 2, line=4, cex=1.3)
    mtext("Intertegular span (mm) ", 1, line=3, cex=1)

    ## 1 persistence
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    legend("topleft", legend="(a)", bty="n", cex=1.2)


    for(i in 1:length(quantiles.B)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.nat.area'] * x +
                        means['phi.B'] * quantiles.B[i] +
                        means['phi.nat.area.B'] * x * quantiles.B[i]),
              from=range(model.input$data$natural)[1],
              to=range(model.input$data$natural)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 2 colonization
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural), las=1,
         ylab="", xlab="")
    legend("topleft", legend="(b)", bty="n", cex=1.2)
    mtext("Colonization", 2, line=4, cex=1.3)

    for(i in 1:length(quantiles.B)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.nat.area'] * x +
                        means['gam.B'] * quantiles.B[i] +
                        means['gam.nat.area.B'] * x * quantiles.B[i]),
              from=range(model.input$data$natural)[1],
              to=range(model.input$data$natural)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## 3 occupancy
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$natural),
         ylab="", xlab="", las=1)
    legend("topleft", legend="(c)", bty="n", cex=1.2)
    mtext("Occupancy", 2, line=4, cex=1.3)

    for(i in 1:length(quantiles.B)){
        curve((inv.logit(means['mu.gam.0'] +
                         means['mu.gam.nat.area'] * x +
                         means['gam.B'] * quantiles.B[i] +
                         means['gam.nat.area.B'] * x *
                         quantiles.B[i]))/
              (1 + inv.logit(means['mu.gam.0'] +
                             means['mu.gam.nat.area'] * x +
                             means['gam.B'] * quantiles.B[i] +
                             means['gam.nat.area.B'] * x *
                             quantiles.B[i]) -
               inv.logit(means['mu.phi.0'] +
                         means['mu.phi.nat.area'] * x +
                         means['phi.B'] * quantiles.B[i] +
                         means['phi.nat.area.B'] * x * quantiles.B[i])),
              from=range(model.input$data$natural)[1],
              to=range(model.input$data$natural)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    mtext("Remnant habitat \n proximity-weighted area", 1, line=5.5, cex=1)

}




