#clear environment
rm(list = ls())

#clear console
cat("\014")

setwd("C:/Chinmay/mika_09_12_2017/data/Mica_whsed_landsat/test")

#import libraries
library(raster)
library(rgdal)
library(animation)


##### read all rasters in a dir and stack them
img <- list.files()
stack <- stack(img)

##### Split raster name to get the date info
split <- sapply(img, function(x) strsplit(x, "_")[[1]], USE.NAMES=FALSE)
date <- split[1,]

shp1<- readOGR("C:/Chinmay/mika_09_12_2017/data/mika_shapefiles/1_wsheds30mwgs84_UTM11N.shp")
shp2<- readOGR("C:/Chinmay/mika_09_12_2017/data/mika_shapefiles/2_wsheds30mwgs84_UTM11N.shp")
shp3<- readOGR("C:/Chinmay/mika_09_12_2017/data/mika_shapefiles/3_wsheds30mwgs84_UTM11N.shp")
shp4<- readOGR("C:/Chinmay/mika_09_12_2017/data/mika_shapefiles/4_wsheds30mwgs84_UTM11N.shp")
shp5<- readOGR("C:/Chinmay/mika_09_12_2017/data/mika_shapefiles/5_wsheds30mwgs84_UTM11N.shp")
shp6<- readOGR("C:/Chinmay/mika_09_12_2017/data/mika_shapefiles/6_wsheds30mwgs84_UTM11N.shp")
shp7<- readOGR("C:/Chinmay/mika_09_12_2017/data/mika_shapefiles/7_wsheds30mwgs84_UTM11N.shp")

#Set delay between frames when replaying
# ani.options(interval=.5)

#ani.options(convert = shQuote(normalizePath('C:/Program Files/ImageMagick-7.0.8-Q16/magick.exe')))

#colors
colfunc <- colorRampPalette(c("brown", "yellow","green", "darkgreen"))
colfunc(10)

fun <- function() {
    plot(shp1, add=T)
    plot(shp2, add=T)
    plot(shp3, add=T)
    plot(shp4, add=T)
    plot(shp5, add=T)
    plot(shp6, add=T)
    plot(shp7, add=T)
}

#visualize in Rstudio
# animate(stack, addfun=fun, zlim=c(0.1, 0.9), main = date, pause=1, n=3, col = colfunc(15))

# Export as GIF file
# saveGIF(animate(stack, addfun=fun, zlim=c(0.1, 0.9), main = date, pause=2, n=3, col = colfunc(15)), movie.name = "animation.gif", ani.width = 800, ani.height = 800, clean = TRUE)

# Export as latex

saveLatex(animate(stack, addfun=fun, zlim=c(0.1, 0.9), main = date, pause=1, n=3, col = colfunc(15)), movie.name = 'animation.gif',ani.width = 800, ani.height = 800,  interval=.8)

