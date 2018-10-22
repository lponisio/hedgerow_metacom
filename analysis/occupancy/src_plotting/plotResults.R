

if(!exists("ms.ms.nimble")){
    load(file=file.path(save.dir,
                        sprintf('runs/%s_%s_%s.Rdata',
                                data.subset,
                                natural.decay, HR.decay)))
}
if(is.list(ms.ms.nimble$samples)){
    samples.4.table <- do.call(rbind, ms.ms.nimble$samples)
} else{
    samples.4.table <- ms.ms.nimble$samples
}

means <- apply(samples.4.table, 2, mean)
se <- apply(samples.4.table, 2, sd)


## pdf.f(f.plotInteractionsHRRemnant.k, file=file.path(save.dir,
##       sprintf("figures/interactions/%s_HRinteractions-k-%s-%s.pdf",
##       data.subset, natural.decay, HR.decay)),
##       width=6, height=11)

pdf.f(plotInteractionsk, file=file.path(save.dir,
       sprintf("figures/interactions/%s_HRinteractions-k-%s-%s.pdf",
       data.subset, natural.decay, HR.decay)),
      width=3, height=9)

pdf.f(f.plotInteractionsFloralDiv, file=file.path(save.dir,
    sprintf("figures/interactions/%s_HRinteractions-fra-%s-%s.pdf",
    data.subset, natural.decay, HR.decay)),
      width=3, height=9)

pdf.f(plotInteractionsB, file=file.path(save.dir,
       sprintf("figures/interactions/%s_HRinteractions-B-%s-%s.pdf",
       data.subset, natural.decay, HR.decay)),
      width=3, height=9)

## ****************************************************************
##  posterior probability table
## ****************************************************************
## posterior probs
if(!exists("ms.ms.model")){
    ms.ms.model <- nimbleModel(code=ms.ms.occ,
                               constants=model.input$constants,
                               data=model.input$data,
                               inits=model.input$inits,
                               check=FALSE,
                               calculate=FALSE)
}

samples.4.table <- samples.4.table[, colnames(samples.4.table) %in%
                    ms.ms.model$getNodeNames(includeData=FALSE,
                                             stochOnly=TRUE)]
## H param > 0
h1 <- apply(samples.4.table,
            2, function(x) sum(x > 0)/length(x))
## H param < 0
h2 <- apply(samples.4.table,
            2, function(x) sum(x < 0)/length(x))
posterior.probs <- cbind(h1,h2)
makeTable()


