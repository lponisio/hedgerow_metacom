## setwd('~/Dropbox/hedgerow_metacom')
setwd('analysis/networks')
source('src/initialize.R')


## load the GIANT model file that tracks all of the latent states
load.dir <- "../../../hedgerow_metacom_saved/occupancy/runs"
load(file.path(load.dir, 'all_2500_350_FALSE.Rdata'))

samples <- do.call(rbind, ms.ms.nimble)
samples.z <- samples[, grepl("Z", colnames(samples))]
samples.z <- colSums(samples.z)/nrow(samples.z)
z.array <- array(samples.z,
                 dim=dim(model.input$data$Z),
                 dimnames=dimnames(model.input$data$Z))
save(z.array,
     file=sprintf('saved/zarray_drop_li_ht%s.Rdata', drop.li.ht))

load(file=sprintf('saved/zarray_drop_li_ht%s.Rdata', drop.li.ht))

specs.year.sp <- vector(mode="list",
                        length=dim(z.array)[2])
specs.year.site <- vector(mode="list",
                          length=dim(z.array)[2])

## connectance of a species  within a year across sites
for(i in 1:dim(z.array)[2]){
    this.year <- z.array[,i,]
    specs.year <- species.lev(this.year)
    ll <- specs.year$`lower level`
    hl <- specs.year$`higher level`
    ll$Year <- dimnames(z.array)[[2]][i]
    ll$Site <- rownames(ll)
    hl$Year <- dimnames(z.array)[[2]][i]
    hl$GenusSpecies <- rownames(hl)
    hl <- merge(hl, traits[,c("GenusSpecies", "r.degree", "MeanITD")],
                key="GenusSpecies")
    ll <- merge(ll, by.site)
    rownames(ll) <- rownames(hl) <- NULL
    specs.year.sp[[i]] <- hl
    specs.year.site[[i]] <- ll
}

specs.years.pol <- do.call(rbind, specs.year.sp)
specs.years.site <- do.call(rbind, specs.year.site)


specs.site.sp <- vector(mode="list",
                        length=dim(z.array)[1])

## connectance of a species  within a site across years
for(i in 1:dim(z.array)[1]){
    this.site <- z.array[i,,]
    specs.site <- species.lev(this.site)
    hl <- specs.site$`higher level`
    hl$Site <- dimnames(z.array)[[1]][i]
    hl$GenusSpecies <- rownames(hl)
    hl <- merge(hl, traits[,c("GenusSpecies", "r.degree", "MeanITD")],
                key="GenusSpecies")
    rownames(hl) <- NULL
    specs.site.sp[[i]] <- hl
}
specs.site.pol <- do.call(rbind, specs.site.sp)


## **********************************************************
## pollinators
## **********************************************************

## anything outputted by specieslevel
ys <- c("tot.int", "weighted.closeness")

xvar.species.form <- paste0("scale(", xvar.species, ")")

formulas.species <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.species.form , collapse="+"),
                           "(1|GenusSpecies)",
                           sep="+")))
})

## within a year, across sites
mod.years.pol <- lapply(formulas.species, function(x){
    lmer(x, data=specs.years.pol)
})

## within a site across years
mod.sites.pol <- lapply(formulas.species, function(x){
    lmer(x, data=specs.site.pol)
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
xvar.site.form <- paste0("scale(", xvar.site, ")")
formulas.site <-lapply(ys, function(x) {
    as.formula(paste(x, "~",
                     paste(paste(xvar.site.form, collapse="+"),
                           "(1|Site)",
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
     file=sprintf('saved/mods/zmets_drop_li_ht%s.Rdata', drop.li.ht))
