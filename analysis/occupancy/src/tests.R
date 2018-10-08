
testOccData <- function(){
    ## ************************************************************
    ## check all the dimensions of the data and the names line up
    ## ************************************************************

    print(paste("natural sites match X matrix",
                all(dimnames(model.input$data$X)$site ==
                                                names(model.input$data$natural))))

    print(paste("floral resource site names matrix matches X matrix",
                all(dimnames(model.input$data$X)$site ==
                                                rownames(model.input$data$fra))))


    print(paste("floral resource year names matrix matches X matrix",
                all(dimnames(model.input$data$X)$year[-10] ==
                                                colnames(model.input$data$fra))))


    print(paste("hedgerow area sites match X matrix",
                all(dimnames(model.input$data$X)$site ==
                                                names(model.input$data$HRarea))))


    print(paste("species degrees match X matrix",
                all(dimnames(model.input$data$X)$species == names(model.input$data$k))))


    print(paste("species body sizes match X matrix",
                all(dimnames(model.input$data$X)$species ==
                                                names(model.input$data$B))))

    ## ************************************************************
    ## check x matrix against raw data
    ## ************************************************************

    raw.check1 <- unique(spec$GenusSpecies[spec$Site == dimnames(model.input$data$X)$site[1]&
                                           spec$Year == dimnames(model.input$data$X)$year[1] &
     spec$Date == min(sr.sched$Date[sr.sched$Site == dimnames(model.input$data$X)$site[1]])])

    model.check1 <- names(model.input$data$X[1,1,1,][ model.input$data$X[1,1,1,] == 1])

    print(paste("random raw, model data match 1", all(raw.check1 %in%
                                                      model.check1)))



    raw.check2 <- unique(spec$GenusSpecies[spec$Site == dimnames(model.input$data$X)$site[5]&
                                         spec$Year == dimnames(model.input$data$X)$year[2] &
  spec$Date == min(sr.sched$Date[sr.sched$Site == dimnames(model.input$data$X)$site[5]])])

    model.check2 <- names(model.input$data$X[5,2,1,][ model.input$data$X[5,2,1,] == 1])

    print(paste("random raw, model data match 2", all(raw.check2 %in%
                                                      model.check2)))


    raw.check3 <- unique(spec$GenusSpecies[spec$Site == dimnames(model.input$data$X)$site[15]&
                                          spec$Year == dimnames(model.input$data$X)$year[5] &
    spec$Date == min(sr.sched$Date[sr.sched$Site == dimnames(model.input$data$X)$site[15]])])

    model.check3 <- names(model.input$data$X[15,5,1,][ model.input$data$X[15,5,1,] == 1])

    print(paste("random raw, model data match 3", all(raw.check3 %in%
                                                      model.check3)))

}
