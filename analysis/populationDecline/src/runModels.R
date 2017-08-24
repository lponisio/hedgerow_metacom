

calcOccArray <-  function(x){
  ## abund
  mean.abund <- mean(x, na.rm=TRUE)
  mean.abund[!is.finite(mean.abund)] <- NA
  cv.abund <- sd(x, na.rm=TRUE)/mean.abund
  cv.abund[!is.finite(cv.abund)] <- NA

  ## occ
  x[x > 1] <- 1
  occ <- mean(x, na.rm=TRUE)
  occ[!is.finite(occ)] <- NA
  cv.occ <- sd(x, na.rm=TRUE)/occ
  cv.occ[!is.finite(cv.occ)] <- NA
  return(c(mean.abund=mean.abund,
           cv.abund=cv.abund,
           occ=occ,
           cv.occ=cv.occ))
}


calcDroughtvar <- function(species.mat,  drought.yrs, traits, spec){
  years <- format(as.Date(dimnames(species.mat)[2]$date), '%Y')
  drought <- apply(species.mat[, years %in% drought.yrs,],
                   c(3,1), calcOccArray)
  nondrought <- apply(species.mat[, !years %in% drought.yrs,],
                      c(3,1), calcOccArray)
  all.dat <- simplify2array(list(drought=drought,
                                 nondrought=nondrought))
  all.dat <- adply(all.dat, 1:(dim(all.dat)[1]))
  colnames(all.dat) <- c("metric",
                         colnames(all.dat)[-c(1,4,5)],
                         "drought", "score")
  all.dat <- merge(all.dat, traits,
                   by.x = "species", by.y = "GenusSpecies")
  all.dat$SiteStatus <- spec$SiteStatus[match(all.dat$site, spec$Site)]
  return(all.dat)
}


## runs models based on a formula, family, response variable and data
runMod <- function(forms,
                   fam,
                   dats,
                   return.sum=FALSE){
  if(fam =="poisson"){
    mod <- glmer(forms,
                 family=fam,
                 data=dats,
                 nAGQ = 10L,
                 control=glmerControl(optimizer="bobyqa",
                   optCtrl=list(maxfun=1e9)))
    ifelse(return.sum,
           return(summary(mod)),
           return(mod))
  }else if(fam =="nbinom"){
    mod <- glmer.nb(forms,
                    data=dats,
                    control=glmerControl(optimizer="bobyqa",
                      optCtrl=list(maxfun=1e9),  tolPwrss=1e-3))
  }else if(fam =="gaussian"){
    mod <- do.call(lmer,
                   list(formula=forms,
                        data=dats))
    ## mod <- lmer(forms,
    ##              data=dats)
  }
  ifelse(return.sum,
         return(summary(mod)),
         return(mod))
}
