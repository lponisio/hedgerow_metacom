plotInteractions <- function(){
    layout(matrix(c(1,1,1, 2:4), nrow=2, byrow=TRUE), heights=c(0.25,1,1,1))
    par(oma=c(6, 7, 2.5, 1),
        mar=c(1, 0, 0, 1.5), cex.axis=1.5)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", bty="n" )
    legend("center", legend=c("min", "5%", "25%", "median", "75%", "95%", "max"),
           pch=16, col=cols, bty="n", cex=1.5, ncol=7)

    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1)
    mtext("Persistence", 2, line=4, cex=1.3)
    mtext("Hedgerow \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.fra'] * quantiles$fra[i] +
                        means['mu.phi.hr.area.fra'] * x * quantiles$fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="a)", bty="n", cex=1.2)

    ## interactions of floral resources and non crop habitat prox
    ## plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural),
    ##      yaxt="n")
    plot(NA, ylim=c(0, 1), xlim=c(-1,1),
         yaxt="n")
    mtext("Non-crop habitat \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.nat.area'] * x +
                        means['mu.phi.fra'] * quantiles$fra[i] +
                        means['mu.phi.nat.area.fra'] * x * quantiles$fra[i]),
              from=-1,
              to=1,
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="b)", bty="n", cex=1.2)

    ## intearction betwen non crop habitat and body size
    plot(NA, ylim=c(0, 1), xlim=c(-1,1), yaxt="n")
    mtext("Non-crop habitat \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.nat.area'] * x +
                        means['mu.phi.B'] * quantiles$B[i] +
                        means['mu.phi.nat.area.B'] * x * quantiles$B[i]),
              from=-1,
              to=1,
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="c)", bty="n", cex=1.2)
}



plotInteractionsLinear <- function(){
    layout(matrix(c(1,1,1, 2:4), nrow=2, byrow=TRUE), heights=c(0.25,1,1,1))
    par(oma=c(6, 7, 2.5, 1),
        mar=c(1, 0, 0, 1.5), cex.axis=1.5)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", bty="n" )
    legend("center", legend=c("min", "5%", "25%", "median", "75%", "95%", "max"),
           pch=16, col=cols, bty="n", cex=1.5, ncol=7)

    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(-5, 5), xlim=range(model.input$data$HRarea), las=1)
    mtext("Persistence", 2, line=4, cex=1.3)
    mtext("Hedgerow \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(means['mu.phi.0'] +
              means['mu.phi.hr.area'] * x +
              means['mu.phi.fra'] * quantiles$fra[i] +
              means['mu.phi.hr.area.fra'] * x * quantiles$fra[i],
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## interactions of floral resources and non crop habitat prox
    plot(NA, ylim=c(-5, 5), xlim=c(-1,1),
         yaxt="n")
    mtext("Non-crop habitat \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(means['mu.phi.0'] +
              means['mu.phi.nat.area'] * x +
              means['mu.phi.fra'] * quantiles$fra[i] +
              means['mu.phi.nat.area.fra'] * x * quantiles$fra[i],
              from=-1,
              to=1,
              col=cols[i],
              lwd=2,
              add=TRUE)
    }

    ## intearction betwen non crop habitat and body size
    plot(NA, ylim=c(-5, 5), xlim=c(-1,1), yaxt="n")
    mtext("Non-crop habitat \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(means['mu.phi.0'] +
              means['mu.phi.nat.area'] * x +
              means['mu.phi.B'] * quantiles$B[i] +
              means['mu.phi.nat.area.B'] * x * quantiles$B[i],
              from=-1,
              to=1,
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
}


plotAllInteractions <- function(){
    layout(matrix(c(1,1,1, 2:7), nrow=3, byrow=TRUE),
           heights=c(0.25,1,1,1,1,1,1))

    par(oma=c(3, 7, 2.5, 1),
        mar=c(5, 0, 0, 1.5), cex.axis=1.5)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", bty="n", ylab="", xlab="")
    legend("center", legend=c("min", "5%", "25%", "median", "75%", "95%", "max"),
           pch=16, col=cols, bty="n", cex=1.5, ncol=7)


    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$k), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    mtext("Floral diet breadth", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * quantiles$HRarea[i] +
                        means['mu.phi.k'] * x +
                        means['mu.phi.hr.area.k'] * x * quantiles$HRarea[i]),
              from=range(model.input$data$k)[1],
              to=range(model.input$data$k)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="Hedgerow proximity (a)", bty="n", cex=1.2, adj=0.5)

    ## intearction betwen non crop habitat and body size
    plot(NA, ylim=c(0, 1), xlim=c(-1,1), yaxt="n", ylab="", xlab="")
    mtext("Body size", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.nat.area'] * x +
                        means['mu.phi.B'] * quantiles$natural[i] +
                        means['mu.phi.nat.area.B'] * x * quantiles$natural[i]),
              from=range(model.input$data$B)[1],
              to=range(model.input$data$B)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="Non-crop habitat proximity (c)", bty="n", cex=1.2)

    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="", yaxt="n")
    ## mtext("Persistence", 2, line=4, cex=1.3)
    mtext("Hedgerow \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.fra'] * quantiles$fra[i] +
                        means['mu.phi.hr.area.fra'] * x * quantiles$fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend=" Floral diversity (d)", bty="n", cex=1.2)


    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$fra), las=1,
         ylab="", xlab="")
    mtext("Floral diversity", 1, line=3, cex=1.3)
    mtext("Colonization", 2, line=4, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * quantiles$HRarea[i] +
                        means['mu.gam.fra'] * x +
                        means['mu.gam.hr.area.fra'] * x * quantiles$HRarea[i]),
              from=range(model.input$data$fra)[1],
              to=range(model.input$data$fra)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="Hedgerow proximity (b)", bty="n", cex=1.2)


    ## interactions of floral resources and non crop habitat prox
    ## plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural),
    ##      yaxt="n")
    plot(NA, ylim=c(0, 1), xlim=c(-1,1),
         yaxt="n", ylab="", xlab="")
    mtext("Non-crop habitat \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.nat.area'] * x +
                        means['mu.phi.fra'] * quantiles$fra[i] +
                        means['mu.phi.nat.area.fra'] * x * quantiles$fra[i]),
              from=-1,
              to=1,
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="Floral diversity (e)", bty="n", cex=1.2)

    ## intearction betwen non crop habitat and body size
    plot(NA, ylim=c(0, 1), xlim=c(-1,1), yaxt="n", ylab="", xlab="")
    mtext("Hedgerow proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.nat.area'] * x +
                        means['mu.gam.B'] * quantiles$B[i] +
                        means['mu.gam.nat.area.B'] * x * quantiles$B[i]),
              from=-1,
              to=1,
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="Body size (f)", bty="n", cex=1.2)

    ## interactions of floral resources and non crop habitat prox
    ## plot(NA, ylim=c(0, 1), xlim=range(model.input$data$natural),
    ##      yaxt="n")
    plot(NA, ylim=c(0, 1), xlim=c(-1,1),
         yaxt="n", ylab="", xlab="")
    mtext("Non-crop habitat \n proximity", 1, line=5, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.nat.area'] * x +
                        means['mu.gam.fra'] * quantiles$fra[i] +
                        means['mu.gam.nat.area.fra'] * x * quantiles$fra[i]),
              from=-1,
              to=1,
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topright", legend="Floral diversity (e)", bty="n", cex=1.2)


}



plotHRInteractions <- function(){
    layout(matrix(c(1,1,1, 2:7), nrow=3, byrow=TRUE),
           heights=c(0.25,1,1,1,1,1,1))

    par(oma=c(3, 7, 2.5, 1),
        mar=c(5, 0, 0, 1.5), cex.axis=1.5)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n", bty="n", ylab="", xlab="")
    legend("center", legend=c("min", "5%", "25%", "median", "75%", "95%", "max"),
           pch=16, col=cols, bty="n", cex=1.5, ncol=7)

    ## 1
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Persistence", 2, line=4, cex=1.3)
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.fra'] * quantiles$fra[i] +
                        means['mu.phi.hr.area.fra'] * x * quantiles$fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Floral diversity (a)", bty="n", cex=1.2)

    ## 2
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="", yaxt="n")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.k'] * quantiles$k[i] +
                        means['mu.phi.hr.area.k'] * x * quantiles$k[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Diet breadth (b)", bty="n", cex=1.2)
    ## 3
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         yaxt="n", ylab="", xlab="")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.phi.0'] +
                        means['mu.phi.hr.area'] * x +
                        means['mu.phi.B'] * quantiles$B[i] +
                        means['mu.phi.hr.area.B'] * x * quantiles$B[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("bottomleft", legend="Body size (c)", bty="n", cex=1.2)

    ## 4
    ## interactions of floral resources and hedgerow proximity
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="")
    mtext("Colonization", 2, line=4, cex=1.3)
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.fra'] * quantiles$fra[i] +
                        means['mu.gam.hr.area.fra'] * x * quantiles$fra[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topleft", legend="Floral diversity (a)", bty="n", cex=1.2)

    ## 5
    plot(NA, ylim=c(0, 1), xlim=range(model.input$data$HRarea), las=1,
         ylab="", xlab="", yaxt="n")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.k'] * quantiles$k[i] +
                        means['mu.gam.hr.area.k'] * x * quantiles$k[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topleft", legend="Diet breadth (b)", bty="n", cex=1.2)

    ## 6
    plot(NA, ylim=c(0, 1),
         xlim=range(model.input$data$HRarea),
         yaxt="n", ylab="", xlab="")
    mtext("Hedgerow proximity", 1, line=3, cex=1.3)

    for(i in 1:length(probs)){
        curve(inv.logit(means['mu.gam.0'] +
                        means['mu.gam.hr.area'] * x +
                        means['mu.gam.B'] * quantiles$B[i] +
                        means['mu.gam.hr.area.B'] * x * quantiles$B[i]),
              from=range(model.input$data$HRarea)[1],
              to=range(model.input$data$HRarea)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }
    legend("topleft", legend="Body size (c)", bty="n", cex=1.2)
}




