
plotComparisons <- function(){
    layout(matrix(1:3, ncol=1))
    par(oma=c(2,1,0.5,1), mar=c(10,2,0,0))
    cols <- rainbow(dim(ms.ms.occ.all$ms.ms$summary)[1])
    for(group in groups){
        if(length(group) > 1){
            for(i in 1:dim(ms.ms.occ.all$ms.ms$summary)[1]){
                this.samp <- ms.ms.occ.all$ms.ms$summary[i,,group]
                xs <- jitter(1:dim(this.samp)[2])
                if(i == 1){
                    plot(x=xs, y=this.samp["mean",], pch=16,
                         col=cols[i],
                         ylim=range(c(ms.ms.occ.all$ms.ms$summary[,'CI95_upp',
                                                                  group],
                                      ms.ms.occ.all$ms.ms$summary[,'CI95_low',
                                                                  group])),
                         xaxt="n",
                         ylab="Estimate",
                         xlab="")
                    abline(h=0,lty=2, col="grey")

                    axis(1, at=1:dim(this.samp)[2],
                         labels=FALSE)
                    text(x=1:dim(this.samp)[2],
                         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
                         labels=group, srt=45, adj=1, xpd=TRUE)
                } else{
                    points(x=xs,
                           y=this.samp["mean",], pch=16, col=cols[i])
                }
                arrows(y1=this.samp['CI95_upp',],
                       y0=this.samp['CI95_low',],
                       x0=xs,
                       code=0, angle=90, length=0.02, lwd=1, col=cols[i])
            }
        }

        legend("topright", legend=dimnames(ms.ms.occ.all$ms.ms$summary)[[1]],
               pch=16, col=cols)
    }
}


plotPosterior <- function(summarys, wanted.order,
                          xlabs, phis, gams){
    plotPhiGam <- function(phi.gam){
        plot(1:length(phi.gam),
             summarys["mean", phi.gam],
             ylim=range(summarys[c("CI95_low", "CI95_upp"), phi.gam]),
             pch=16, las=1,
             ylab="", xlab="", xaxt="n",
             xlim= range(c(1,length(phi.gam)) +  c(-0.25, 0.25)))
        arrows(y1=summarys['CI95_upp', phi.gam],
               y0=summarys['CI95_low', phi.gam],
               x0=1:length(wanted.order),
               code=0, angle=90, length=0.02, lwd=1)
        abline(h=0, col="gray45", lty=2)
    }
    ## layout(matrix(c(1,2,3,4, rep(5, 4), rep(6, 4)), nrow=3,
    ## byrow=TRUE))
    layout(matrix(1:2, nrow=2))
    par(oma=c(0, 4, 1, 1),
        mar=c(2, 4, 5.5, 1), cex.axis=1.2)
    ## phis <- paste("mu.phi", wanted.order,
    ##               sep=".")
    ## gams <- paste("mu.gam", wanted.order,
    ##               sep=".")

    plotPhiGam(phis)
    legend("topright", legend="(a)", bty="n", cex=1.2)
    mtext(text="Persistence", 2,
          line=3.75, cex=1.3)
    old.mar <- par("mar")
    par(mar = old.mar + c(5.5,0,-5.5,0))

    plotPhiGam(gams)
    legend("topright", legend="(b)", bty="n", cex=1.2)
    axis(1, at=1:length(phis),
         labels=FALSE)
    text(x=1:length(phis),
         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
         labels=xlabs, srt=45, xpd=TRUE, adj=1, cex=0.95)
    mtext(text="Colonization", 2,
          line=3.75, cex=1.3)
    mtext(text="Posterior model estimate", 2,
          line=6, cex=1.3, at=max(summarys['CI95_upp', gams]) + 0.5)
}


plotVariables <- function(){
    layout(matrix(1:2, ncol=1))
    par(oma=c(2, 2, 1, 1),
        mar=c(4, 4, 2, 1))
    hist(model.input$data$HRarea, main="", xlab="", ylab="", las=1)
    abline(v=0, lty=2, col="red")
    legend("topright", legend="(a)", bty="n", cex=1.2)
    mtext("Hedgerow area/proximity", 1, line=3, cex=0.9)
    mtext(text="Frequency", 2,
          line=3.75, cex=0.9)
    hist(model.input$data$natural, main="", xlab="", ylab="", las=1)
    abline(v=0, lty=2, col="red")
    legend("topright", legend="(b)", bty="n", cex=1.2)
    mtext("Remnant habitat area/proximity", 1, line=3, cex=0.9)
    mtext(text="Frequency", 2,
          line=3.75, cex=0.9)

}

plotVariablesWrapper <- function(){
    checkDirExists(file.path(save.dir, 'figures/variables'))
    pdf.f(plotVariables,
          file=file.path(save.dir, 'figures/variables',
                         sprintf('%s%s.pdf', natural.decay, HR.decay)),
          height= 6, width=3)

}



## plotVariables <- function(){
##     layout(matrix(1:2, ncol=1))
##     par(oma=c(2, 2, 1, 1),
##         mar=c(4, 4, 2, 1))

    ## plot(density(model.input$data$HRarea[names(model.input$data$HRarea) %in% controls]),
    ##      main="", xlab="", ylab="",
    ##      xlim=range(model.input$data$HRarea) + c(0,1), las=1, col="gray45")
    ## points(density(model.input$data$HRarea[names(model.input$data$HRarea) %in% hedgerows]),
    ##        type="l")
    ## hist(model.input$data$HRarea, main="", xlab="", ylab="", las=1)
    ## abline(v=0, lty=2, col="red")
    ## legend("topright", legend="(a)", bty="n", cex=1.2)
    ## mtext("Hedgerow area/proximity", 1, line=3, cex=0.9)
    ## mtext(text="Frequency", 2,
    ##       line=3.75, cex=0.9)
    ## legend("topright", legend=c("Hedgerows", "Field margins"),
    ##        col=c("black", "gray45"), lty=1, bty="n", cex=0.8)

    ## plot(density(model.input$data$natural[names(model.input$data$natural) %in% controls]),
    ##      main="", xlab="", ylab="", las=1,
    ##      xlim=range(model.input$data$natural) + c(0,1), col="gray45")
    ## points(density(model.input$data$natural[names(model.input$data$natural) %in% hedgerows]),
    ##        type="l", lty=1)
    ## hist(model.input$data$natural, main="", xlab="", ylab="", las=1)
    ## abline(v=0, lty=2, col="red")
    ## legend("topright", legend="(b)", bty="n", cex=1.2)
    ## mtext("Remnant habitat area/proximity", 1, line=3, cex=0.9)
    ## mtext(text="Frequency", 2,
    ##       line=3.75, cex=0.9)

    ## plot(density(by.site$Div[by.site$Site %in% hedgerows]),
    ##      main="", xlab="", ylab="", las=1)
    ## points(density(by.site$Div[by.site$Site %in% controls]),
    ##        type="l", col="gray45")
    ## abline(v=mean(by.site$Div, na.rm=TRUE), lty=2)
    ## hist(model.input$data$fra, main="", xlab="", ylab="", las=1)
    ## abline(v=mean(model.input$data$fra, na.rm=TRUE), lty=2, col="red")
    ## legend("topright", legend="(c)", bty="n", cex=1.2)
    ## mtext("Floral diversity", 1, line=3, cex=0.9)
    ## mtext(text="Frequency", 2,
    ##       line=3.75, cex=0.9)

    ## hist(all.traits$BodyLength, main="", xlab="", las=1, ylab="")
    ## abline(v=mean(all.traits$BodyLength, na.rm=TRUE), lty=2)

    ## hist(model.input$data$B, main="", xlab="", las=1, ylab="")
    ## abline(v=mean(model.input$data$B, na.rm=TRUE), lty=2, col="red")
    ## legend("topright", legend="(d)", bty="n", cex=1.2)
    ## mtext("Body size (mm)", 1, line=3, cex=0.9)
    ## mtext(text="Frequency", 2,
    ##       line=3.75, cex=0.9)

## }

