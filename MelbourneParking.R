# adjust path
setwd("~/NetworkSplines")

# load required packages
library(spatstat)
library(igraph)
library(Matrix)
library(MASS)
library(VCA)
library(splines)
library(dplyr)
library(lubridate)
library(readxl)

# load functions
source(file = "Functions.R")


#### geometric network ####

# edges and coordinates of vertices
E <- read_xlsx("Data/EdgesMelbourne.xlsx")
V <- read_xlsx("Data/VerticesMelbourne.xlsx") %>% mutate(Lat = as.numeric(Lat), Lon = as.numeric(Lon))

# change coordinate system
min.lon <- min(V$Lon)
min.lat <- min(V$Lat)
V$Lon <- V$Lon - min.lon
V$Lat <- V$Lat - min.lat

# create linnet object and rescale with units meters
P <- ppp(x = V$Lon, y = V$Lat, 
        window = owin(xrange = c(min(V$Lon), max(V$Lon)), yrange = c(min(V$Lat), max(V$Lat))))
L <- linnet(vertices = P, edges = as.matrix(E[, 2:3]))
s <- L$dpath[12, 94]/1900
L <- spatstat::rescale(L, s = s, unitname = "Meters")

# choose parameters and augment L
delta <- 50
h <- 5
r <- 2
L <- augment.linnet(L, delta, h)

# read parking data
data <- readRDS("Data/data_2019_clean.rds")
lots <- readRDS("Data/Parking_Lots.rds") %>% filter(sensor == 1) %>%
  select(lon, lat, StreetMarker) %>%
  mutate(lon = (lon - min.lon)/s, lat = (lat - min.lat)/s)
data <- left_join(data, lots, by = "StreetMarker")

L.ppp <- ppp(x = lots$lon , y = lots$lat, L$window, marks = lots$StreetMarker)
L.psp <- as.psp(L)
projection <- project2segment(L.ppp, L.psp)

# choose parking lots which are less than 20 meters away from L
seg <- projection$mapXY[which(projection$d < 20)]
tp <- projection$tp[which(projection$d < 20)]
marks <- projection$Xproj$marks[which(projection$d < 20)]

# create L.lpp object with nearby parking lots
L.lpp <- as.lpp(seg = seg, tp = tp, L = L, marks = marks)

# intensity estimate of parking lots
intens.lots <- intensity.pspline.lpp(L.lpp, r)

# only where intensity is at least 0.1
intens.lots.adj <- intens.lots
intens.lots.adj$v[which(intens.lots$v < 0.1, arr.ind = TRUE)] = NA

# plot parking lots on the geometric network
pdf(file = "Plots/MelbourneLots.pdf", width = 10, height = 8)
par(mar=c(0, 1, 2, 1), cex = 2.4)
plot(L.lpp, use.marks = FALSE, pch = 16, cols = "red", lwd = 3, legend = FALSE,
     main = "On-Street Parking Lots in the \n Western CBD of Melbourne, Australia")
dev.off()

# plot intensity of parking lots
min.intens <- min(intens.lots$v[-which(is.na(intens.lots$v))])
max.intens <- max(intens.lots$v[-which(is.na(intens.lots$v))])
pdf(file = "Plots/MelbourneIntensityLots.pdf", width = 10, height = 8)
par(mar=c(0, 1, 2, 1), cex = 2.4)
plot(intens.lots, box = TRUE, 
     main = "Intensity of On-Street Parking Lots in \n the Western CBD of Melbourne, Australia")
dev.off()



# parking data ----

# only consider parking lots which are located on the map
data.CBD <- filter(data, is.element(data$StreetMarker, L.lpp$data$marks), h.start >= 8, h.start < 20)
marks <- select(data.CBD, h.start, weekday)

# observed processes
L.lpp.parking <- as.lpp(x = data.CBD$lon, y = data.CBD$lat, L = L, marks = marks)

#
intens.morning <- intensity.pspline.lpp(L.lpp.morning, r)/52

