## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')
source('src/initialize.R')

## **********************************************************
## pollinators
## **********************************************************

## anything outputted by specieslevel
ys <- c("k", "weighted.betweenness")
formulas.species <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.species, collapse="+"),
                           "(1|Site)",
                           "(1|GenusSpecies)",
                           sep="+")))
})

## within a site across years
mod.sites.pol <- lapply(formulas.species, function(x){
    lmer(x, data=specs.site.pol)
})

## within a year, across sites
mod.years.pol <- lapply(formulas.species, function(x){
    lmer(x, data=specs.years.pol)
})


names(mod.years.pol) <- names(mod.sites.pol) <- ys

## floral degree positivly related to all measures of network
## importance in space/time
print("metacommunity pollinator spatial network")
print(lapply(mod.years.pol, summary))
print("metacommunity pollinator temporal network")
print(lapply(mod.sites.pol, summary))

## **********************************************************
## sites
## **********************************************************

formulas.site <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.site, collapse="+"),
                           "(1|Site)",
                           "(1|GenusSpecies)",
                           sep="+")))
})

## within a year across sites
mod.years.site <- lapply(formulas.site, function(x){
    lmer(x, data=specs.years.site)
})

## name them the same as the pollinators
names(mod.years.site)  <- ys

print("metacommunity patch network")
print(lapply(mod.years.site, summary))

save(mod.years.pol, mod.sites.pol, mod.years.site,
     file=sprintf('saved/mods/specmetrics_drop_li_ht%s.Rdata', drop.li.ht))



load('~/Dropbox/hedgerow_metacom_saved/occupancy/runs/Z_all_2500_350.Rdata')
samples <- do.call(rbind, ms.ms.nimble)

samples.z <- samples[, grepl("Z", colnames(samples))]
samples.z <- colSums(samples.z)/nrow(samples.z)


dim1 <- sapply(strsplit(names(samples.z), ","), function(x) x[1])
dim1 <- as.numeric(gsub("Z\\[", "", dim1))
dim2 <- as.numeric(sapply(strsplit(names(samples.z), ","), function(x) x[2]))
dim3 <- sapply(strsplit(names(samples.z), ","), function(x) x[3])
dim3 <- as.numeric(gsub("\\]", "", dim3))

z.array <- array(samples.z,
                 dim=dim(model.input$data$Z),
                  dimnames=dimnames(model.input$data$Z))
