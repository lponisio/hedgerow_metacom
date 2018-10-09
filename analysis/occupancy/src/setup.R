
## ************************************************************
## specify the parameters to be monitored
## ************************************************************

getParams <- function(){
    c('mu.p.0',
      'sigma.p.0',
      'mu.p.day.1',
      'sigma.p.day.1',
      'mu.p.day.2',
      'sigma.p.day.2',

      'mu.phi.0',
      'sigma.phi.0',
      'mu.phi.hr.area',
      'sigma.phi.hr.area',
      'mu.phi.nat.area',
      'sigma.phi.nat.area',
      'mu.phi.fra',
      'sigma.phi.fra',
      'phi.k',

      'phi.B',
      'phi.nat.area.fra',
      'phi.hr.area.fra',
      'phi.nat.area.k',
      'phi.hr.area.k',
      'phi.nat.area.B',
      'phi.hr.area.B',

      'mu.gam.0',
      'sigma.gam.0',
      'mu.gam.hr.area',
      'sigma.gam.hr.area' ,
      'mu.gam.nat.area',
      'sigma.gam.nat.area',
      'mu.gam.fra',
      'sigma.gam.fra',

      'gam.k',
      'gam.B',
      'gam.hr.area.fra',
      'gam.nat.area.fra',
      'gam.hr.area.k',
      'gam.nat.area.k',
      'gam.hr.area.B',
      'gam.nat.area.B'


      ## ## site level effects
      ## 'phi.nat.area',
      ## 'phi.hr.area',
      ## 'phi.hr.area.fra',
      ## 'phi.nat.area.fra',
      ## 'phi.hr.area.k'

      )
}

getInits <- function(nsp){
    list(mu.p.0 = rnorm(1),
         sigma.p.0 = runif(1),
         mu.p.day.1 = rnorm(1),
         sigma.p.day.1 = runif(1),
         mu.p.day.2 = rnorm(1),
         sigma.p.day.2= runif(1),

         mu.phi.0 = rnorm(1),
         sigma.phi.0 = runif(1),

         mu.phi.hr.area = rnorm(1),
         sigma.phi.hr.area = runif(1),
         mu.phi.nat.area = rnorm(1),
         sigma.phi.nat.area = runif(1),
         mu.phi.fra = rnorm(1),
         sigma.phi.fra = runif(1),

         mu.gam.0 = rnorm(1),
         sigma.gam.0 = runif(1),

         mu.gam.hr.area = rnorm(1),
         sigma.gam.hr.area = runif(1),
         mu.gam.nat.area = rnorm(1),
         sigma.gam.nat.area = runif(1),
         mu.gam.fra = rnorm(1),
         sigma.gam.fra = runif(1),

         p.0 = rnorm(nsp),
         p.day.1 = rnorm(nsp),
         p.day.2 = rnorm(nsp),

         phi.0 = rnorm(nsp),
         phi.hr.area = rnorm(nsp),
         phi.nat.area = rnorm(nsp),
         phi.fra = rnorm(nsp),

         phi.k = rnorm(1),
         phi.B = rnorm(1),
         phi.hr.area.fra = rnorm(1),
         phi.nat.area.fra = rnorm(1),
         phi.nat.area.k = rnorm(1),
         phi.hr.area.k = rnorm(1),
         phi.nat.area.B = rnorm(1),
         phi.hr.area.B = rnorm(1),

         gam.0 = rnorm(nsp),
         gam.hr.area = rnorm(nsp),
         gam.nat.area = rnorm(nsp),
         gam.fra = rnorm(nsp),

         gam.k = rnorm(1),
         gam.B = rnorm(1),
         gam.hr.area.fra = rnorm(1),
         gam.nat.area.fra = rnorm(1),
         gam.hr.area.k = rnorm(1),
         gam.nat.area.k = rnorm(1),
         gam.hr.area.B = rnorm(1),
         gam.nat.area.B = rnorm(1)
         )

}



## random objects needed for plotting
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
