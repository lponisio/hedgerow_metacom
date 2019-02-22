
plotLandscapeCovers <- function(){
    ## prep map and site data
    all.sites.pt <- spTransform(all.sites.pt,
                                CRS("+init=epsg:4326"))

    landcover.nat <- spTransform(landcover.nat,
                                 CRS("+init=epsg:4326"))

    bbox.sites <- bbox(all.sites.pt)
    bbox.nat <- bbox(landcover.nat)

    ## no idea why it doesn't work the first time
    bbox.all <- matrix(c(bbox.sites[1,1],
                         bbox.sites[2,1],
                         bbox.sites[1,2],
                         bbox.sites[2,2]),
                       ncol=2)

    sys <- try(gmap(x=bbox.all,
                    scale=2,
                    type="satellite"))

    sys <- gmap(x=bbox.all,
                scale=2,
                type="satellite", zoom=10)

    all.sites.pt <- spTransform(all.sites.pt,
                                sys@crs)

    landcover.nat <- spTransform(landcover.nat,
                                 sys@crs)


    new.hr <- all.sites.pt[grepl("new", all.sites.pt@data$df0),]
    surveyed.hr <- all.sites.pt[all.sites.pt@data$df0 %in%
                                unique(spec$Site[spec$SiteStatus %in%
                                                 c("maturing", "mature")]),]
    surveyed.control <- all.sites.pt[all.sites.pt@data$df0 %in%
                                     unique(spec$Site[spec$SiteStatus =="control"]),]

    dims <- bbox(sys)
    plot(sys)
    rect(xleft=dims[1,1],ybottom=dims[2,1],
         xright=dims[1,2],ytop=dims[2,2],
         col= rgb(1,1,1, alpha=0.3))

    natural.cover <- crop(landcover.nat,  dims)
    only.natural.cover <-  natural.cover[!natural.cover@data$VegName
                                         %in% kinda.natural,]
    non.native.natural.cover <- natural.cover[natural.cover@data$VegName
                                              %in% kinda.natural,]

    colors <-
        add.alpha(rainbow((length(unique(natural.cover@data$VegName))
            +1 )), 0.5)
    names(colors) <- unique(natural.cover@data$VegName)

    for(veg.type in unique(natural.cover@data$VegName)){
        plot(natural.cover[natural.cover@data$VegName == veg.type,],
             col=colors[veg.type], add=TRUE)
    }

    points(new.hr, col="darkturquoise", pch=15)
    points(surveyed.hr, col="dodgerblue", pch=16)
    points(surveyed.control, col="navy", pch=17)

    points(new.hr, col="black", pch=0)
    points(surveyed.hr, col="black", pch=1)
    points(surveyed.control, col="black", pch=2)

    leg.names <- c("Alkali Bulrush",
                   "Blackberry",
                   "Blue Oak Alliance",
                   "Bullrush/Cattail Wetland",
                   "Bulrush/Cattail Fresh Water Marsh",
                   "Annual Grasslands",
                   "Juniper",
                   "Canyon Live Oak",
                   "Wet Meadow Grasses",
                   "Chamise, Ceanothus",
                   "Chamise",
                   "Coyote Brush",
                   "Wetland Grasses",
                   "Eucalyptus",
                   "Evergreen Shrubland",
                   "Foothill Pine/Chaparral",
                   "Foothill Pine",
                   "Fremont Cottonwood/Valley Oak/Willow",
                   "Valley Oak Riparian Association",
                   "Live Oak/Blue Oak",
                   "Live Oak",
                   "Deciduous Shrubland",
                   "Lotus scoparius",
                   "Mixed Fremont Cottonwood",
                   "Mixed Oak",
                   "Mixed Willow",
                   "Perennial pepperweed",
                   "Saltgrass",
                   "Scrub Oak Chaparral",
                   "Sparse Bush Lupine/Annual Grasses",
                   "Sparse Juniper/Live Oak/Bay/Buckeye",
                   "Tamarisk",
                   "Toyon/Annual Grasses Savanna",
                   "Upland Annual Grasslands & Forbs",
                   "Valley Oak/Cottonwood/Riparian Forest",
                   "Valley Oak",
                   "Valley Oak - Riparian",
                   "Vernal Pool Complex",
                   "White Alder/Riparian Forest",
                   "White Leaf Manzanita/Xeric Serpentine")

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n")
    legend("left", col=colors, pch=16,
           legend=leg.names, ncol=3, bty="n", cex=0.5)
    legend("left", col="black", pch=1,
           legend=leg.names, ncol=3, bty="n",cex=0.5)

    only.2.cols <- add.alpha(c("red", "violetred"), 0.5)
    plot(sys)
    plot(non.native.natural.cover, col=only.2.cols[1], add=TRUE)
    plot(only.natural.cover, col=only.2.cols[2], add=TRUE)

    points(new.hr, col="darkturquoise", pch=15)
    points(surveyed.hr, col="dodgerblue", pch=16)
    points(surveyed.control, col="navy", pch=17)

    points(new.hr, col="black", pch=0)
    points(surveyed.hr, col="black", pch=1)
    points(surveyed.control, col="black", pch=2)

    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n")
    legend("left", col=only.2.cols, pch=16,
           legend=c("Predominantly non-native vegetation",
                    "Predominantly native vegetation"),
           ncol=1, bty="n")
    legend("left", col="black", pch=1,
           legend=c("Predominantly non-native vegetation",
                    "Predominantly native vegetation"),
           ncol=1, bty="n")


    plot(NA, ylim=c(0,1), xlim=c(0,1), xaxt="n", yaxt="n")
    legend("left",
           col=c("darkturquoise","dodgerblue", "navy"),
           pch=c(15, 16, 17),
           legend=c("Unsurveyed hedgerow",
                    "Surveyed hedgerow",
                    "Surveyed field margin"),
           ncol=3, bty="n")
    legend("left", col="black",
           pch=c(0, 1, 2),
           legend=c("Unsurveyed hedgerow",
                    "Surveyed hedgerow",
                    "Surveyed field margin"),
           ncol=3, bty="n")

}


plotLandscapeCoversWithBuffers <- function(){
    ## prep map and site data
    all.sites.pt <- spTransform(all.sites.pt,
                                CRS("+init=epsg:4326"))

    landcover.nat <- spTransform(landcover.nat,
                                 CRS("+init=epsg:4326"))

    bbox.sites <- bbox(all.sites.pt)
    bbox.nat <- bbox(landcover.nat)

    ## no idea why it doesn't work the first time
    bbox.all <- matrix(c(bbox.sites[1,1],
                         bbox.sites[2,1],
                         bbox.sites[1,2],
                         bbox.sites[2,2]),
                       ncol=2)
    sys <- try(gmap(x=bbox.all,
                    scale=2,
                    type="satellite"))
    sys <- try(gmap(x=bbox.all,
                    scale=2,
                    type="satellite", zoom=10))
    all.sites.pt <- spTransform(all.sites.pt,
                                sys@crs)
    landcover.nat <- spTransform(landcover.nat,
                                 sys@crs)
    new.hr <- all.sites.pt[grepl("new", all.sites.pt@data$df0),]
    surveyed.hr <- all.sites.pt[all.sites.pt@data$df0 %in%
                                unique(spec$Site[spec$SiteStatus %in%
                                                 c("maturing", "mature")]),]
    surveyed.control <- all.sites.pt[all.sites.pt@data$df0 %in%
                                     unique(spec$Site[spec$SiteStatus =="control"]),]

    dims <- bbox(sys)
    plot(sys)
    rect(xleft=dims[1,1],ybottom=dims[2,1],
         xright=dims[1,2],ytop=dims[2,2],
         col= rgb(1,1,1, alpha=0.3))
    natural.cover <- crop(landcover.nat,  dims)
    poly.cols <- add.alpha("blue",
                           alpha=0.2)
    plot(natural.cover,
         col=poly.cols, add=TRUE)

    buffers <- gBuffer(rbind(surveyed.hr, surveyed.control),
                       width = 1000)
    plot(buffers, col = add.alpha("yellow", .35), add = TRUE)
    buffers2 <- gBuffer(rbind(surveyed.hr, surveyed.control),
                        width = 7000)
    plot(buffers2, col = add.alpha("red", .35), add = TRUE)

    points(new.hr, col="darkturquoise", pch=15, cex=0.5)
    points(surveyed.hr, col="dodgerblue", pch=16, cex=0.5)
    points(surveyed.control, col="navy", pch=17, cex=0.5)
    points(new.hr, col="black", pch=0, cex=0.5)
    points(surveyed.hr, col="black", pch=1, cex=0.5)
    points(surveyed.control, col="black", pch=2, cex=0.5)
}
