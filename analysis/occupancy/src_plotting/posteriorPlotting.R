
plotPosteriors <- function(){
    if(is.list(ms.ms.nimble$samples)){
        all.samples <- do.call(rbind, ms.ms.nimble$samples)
    } else{
        all.samples <- ms.ms.nimble$samples
    }

    nimble.summary <- apply(all.samples, 2, function(x){
        means <- mean(x)
        CI95_upp <- 1.96*sd(x) + means
        CI95_low <-  means - 1.96*sd(x)
        return(c(mean=means, CI95_upp=CI95_upp, CI95_low=CI95_low ))
    })
    mus <- nimble.summary[, to.plot]

    if(include.int == "no_noncrop"){
        to.plot <- to.plot[!grepl("nat", to.plot)]
        xlabs <- xlabs[!grepl("Non-crop", xlabs)]
    }

    f <- function() {plotPosterior(mus, to.plot, xlabs, phis, gams)}
    checkDirExists(file.path(save.dir, "figures/posterior"))
    pdf.f(f,
          file=file.path(save.dir,
                         sprintf('figures/posterior/%s_mus_bees_%s_%s_%s.pdf',
                                 data.subset,
                                 natural.decay, HR.decay,
                                 include.int)),
          height=7, width=12)

    ## save(mus, file=file.path(save.dir,
    ##                          sprintf('runs/%s_mus_bees_%s_%s_%s.Rdata',
    ##                                  data.subset,
    ##                                  natural.decay, HR.decay,
    ##                                  include.int)))


    all.samples <- all.samples[, colnames(all.samples) %in%
                                 ms.ms.model$getNodeNames(includeData=FALSE,
                                                          stochOnly=TRUE)]
    return(all.samples)
}

makeTable <- function(){
    probs.4.table.phis <- round(posterior.probs[rownames(posterior.probs) %in%
                                                phis,],3)[phis,]
    probs.4.table.gams <- round(posterior.probs[rownames(posterior.probs) %in%
                                                gams,], 3)[gams,]
    means.4.table.phis <- round(means[names(means) %in% phis],3)[phis]
    means.4.table.gams <- round(means[names(means) %in% gams],3)[gams]
    ses.4.table.phis <- round(se[names(se) %in% phis],3)[phis]
    ses.4.table.gams <- round(se[names(se) %in% gams],3)[gams]
    names.4.table <- gsub("\\n", "", xlabs)

    mat.4.table <- cbind(paste0(means.4.table.phis, "(", ses.4.table.phis, ")"),
                         probs.4.table.phis,
                         paste0(means.4.table.gams, "(", ses.4.table.gams, ")"),
                         probs.4.table.gams)

    rownames(mat.4.table) <- names.4.table
    checkDirExists(file.path(save.dir, "table"))
    write.table(mat.4.table,
                sep=" & ",
                file= file.path(save.dir,
                                sprintf('table/%s_post_probs_nimble_bees_%s_%s.txt',
                                        data.subset,
                                        natural.decay, HR.decay)))
}

checkChains <- function(){
    checkDirExists(file.path(save.dir, "figures/chains"))
    runMCMCcheckChains(ms.ms.nimble$samples,
                       f.path= file.path(save.dir,
                                         'figures/chains'),
                       natural.decay, include.int,
                       data.subset, num.samps=10000)
}



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
