
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
                          xlabs){
    plotPhiGam <- function(phi.gam){
        plot(1:length(wanted.order),
             summarys["mean", phi.gam],
             ylim=range(summarys[c("CI95_low", "CI95_upp"), phi.gam]),
             pch=16,
             ylab="", xlab="", xaxt="n",
             xlim= range(c(1,length(wanted.order)) +  c(-0.5, 0.5)))
        arrows(y1=summarys['CI95_upp', phi.gam],
               y0=summarys['CI95_low', phi.gam],
               x0=1:length(wanted.order),
               code=0, angle=90, length=0.02, lwd=1)
        abline(h=0, col="grey", lty=2)
    }
    layout(matrix(1:2, ncol=1))
    par(oma=c(0, 3, 2, 1),
        mar=c(2, 4, 5.5, 1), cex.axis=1.5)
    phis <- paste("mu.phi", wanted.order,
                  sep=".")
    gams <- paste("mu.gam", wanted.order,
                  sep=".")
    plotPhiGam(phis)
    mtext(text="Persistence", 2,
          line=3, cex=1.5)
    old.mar <- par("mar")
    par(mar = old.mar + c(5.5,0,-5.5,0))
    plotPhiGam(gams)
    axis(1, at=1:length(wanted.order),
         labels=FALSE)
    text(x=1:length(wanted.order),
         y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
         labels=xlabs, srt=45, xpd=TRUE, adj=1, cex=1)
    mtext(text="Colonization", 2,
          line=3, cex=1.5)
    mtext(text="Posterior model estimate", 2,
          line=5, cex=1.5, at=max(summarys['CI95_upp', gams]) + 0.5)
}

