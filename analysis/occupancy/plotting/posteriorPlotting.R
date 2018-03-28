
wanted.order <- c("hr.area",
                  "nat.area",
                  "fra",
                  "k",
                  "B",
                  "hr.area.fra",
                  "nat.area.fra",
                  "hr.area.k",
                  "nat.area.k",
                  "hr.area.B",
                  "nat.area.B")


phis <- paste("phi", wanted.order,
              sep=".")
phis <- paste(c(rep("mu.", 3), rep("", 8)), phis, sep="")
gams <- paste("gam", wanted.order,
              sep=".")
gams <- paste(c(rep("mu.", 3), rep("", 8)), gams, sep="")

to.plot <- c(phis, gams)

xlabs <- c("Hedgerow \n area/proximity",
           "Non-crop habitat \n area/proximity",
           "Floral diversity",
           "Floral diet breadth",
           "Body size",
           "Hedgerow \n area/proximity*\n floral diversity",
           "Non-crop \n area/proximity*\n floral diversity",
           "Hedgerow \n area/proximity*\n floral diet breadth",
           "Non-crop \n area/proximity*\n floral diet breadth",
           "Hedgerow \n area/proximity*\n body size",
           "Non-crop \n area/proximity*\n body size")

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
    ## mus <- nimble.summary[,grep("^mu", colnames(nimble.summary))]
    mus <- nimble.summary[, to.plot]

    if(include.int == "no_noncrop"){
        to.plot <- to.plot[!grepl("nat", to.plot)]
        xlabs <- xlabs[!grepl("Non-crop", xlabs)]
    }

    f <- function() {plotPosterior(mus, to.plot, xlabs, phis, gams)}

    pdf.f(f,
          file=file.path(save.dir,
                         sprintf('figures/posterior/%s_mus_bees_%s_%s.pdf',
                                 data.subset, natural.decay, include.int)),
          height=7, width=12)

    save(mus, file=file.path(save.dir,
                             sprintf('runs/%s_mus_bees_%s_%s.Rdata',
                                     data.subset, natural.decay, include.int)))


    all.samples <- all.samples[, colnames(all.samples) %in%
                                 ms.ms.model$getNodeNames(includeData=FALSE,
                                                          stochOnly=TRUE)]
}

makeTable <- function(){
    ## phis <- paste("phi", wanted.order,
    ##               sep=".")
    ## gams <- paste("gam", wanted.order,
    ##               sep=".")

    probs.4.table.phis <- round(posterior.probs[rownames(posterior.probs) %in%
                                                phis,],3)[phis,]
    probs.4.table.gams <- round(posterior.probs[rownames(posterior.probs) %in%
                                                gams,], 3)[gams,]

    rownames(probs.4.table.phis) <- rownames(probs.4.table.gams) <- xlabs

    write.table(cbind(probs.4.table.phis[,-2], probs.4.table.gams[,-2]),
                sep=" & ",

                file= file.path(save.dir,
                     sprintf('table/%s_post_probs_nimble_bees_%s_%s.txt',
                                data.subset, natural.decay, include.int)))
}

checkChains <- function(){
    runMCMCcheckChains(ms.ms.nimble$samples,
                       f.path= file.path(save.dir,
                                         'figures/chains'),
                       natural.decay, include.int,
                       data.subset)
}

