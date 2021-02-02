library(maps)
library(mapdata)
library(mapproj) # install.packages("ggplot2")
library(maptools)  #for shapefiles
library(scales)  #for transparency
locs <- read.csv("coordinates.csv")   #my data for sampling sites, contains a column of "lat" and a column of "lon" with GPS points in decimal degrees

# setEPS()
# postscript("mapUT.eps")
# map("state", "Utah")
# dev.off()

# Plot and export to EPS
setEPS()
postscript("../map.eps")
# map("worldHires","usa", xlim=c(-124.7, -67.1), ylim = c(25.2, 49.4), col="lemonchiffon", fill=TRUE)
# map.cities(x = data(us.cities), label = NULL, minpop=200000, maxpop=Inf,capitals=0, pch=19, cex=0.5,add=TRUE)
map("state", col="gray98", fill=TRUE) #, projection="mercator"
#plot my cities
points(locs$lon, locs$lat, pch=19, col="red", cex=0.8)
dev.off()
