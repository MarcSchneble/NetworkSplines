# adjust path
setwd("~/NetworkSplines")
rm(list = ls())

# load required packages
library(spatstat)
library(igraph)
library(Matrix)
library(MASS)
library(splines)
library(dplyr)
library(lubridate)
library(readxl)
library(tidyr)
library(ggmap)

# load functions
source(file = "Functions.R")
register_google("AIzaSyDX0CVJsIDBJF8NVFPAH84oLWPfvPa335Y")

data <- read.csv(file = "Data/CarCrashesMontgomery.csv") %>% as_tibble() %>%
  mutate(Latitude = as.numeric(substring(sub("\\,.*", "", Location), first = 2)),
         Longitude = as.numeric(gsub("[)]" ,"" , sub("^\\S+", "", Location))))

data.highways <- filter(data, Route.Type %in% c("Maryland (State)", "Interstate (State)", "US (State)"))

#g <- ggmap(get_map(location = c(-77.1, 39.), zoom = 13, maptype = "roadmap")) +
#  geom_point(data = data.highways, aes(x = Longitude, y = Latitude))

V <- read_excel("Data/VerticesMontgomery.xlsx") %>% filter(!is.na(lat), is.na(remove)) %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon))
E <- read_excel("Data/EdgesMontgomery.xlsx") %>% filter(!is.na(from))

E$from.lat <- V$lat[match(E$from, V$Id)]
E$from.lon <- V$lon[match(E$from, V$Id)]
E$to.lat <- V$lat[match(E$to, V$Id)]
E$to.lon <- V$lon[match(E$to, V$Id)]

left <- -77.2 
bottom <- 38.94 
right <- -77.05 
top <- 39.04

g <- ggmap(get_map(location = c(left = left, bottom = bottom, right = right, top = top), maptype = "roadmap", scale = 2)) +
  geom_segment(data = E, aes(x = from.lon, y = from.lat, xend = to.lon, yend = to.lat), size = 1.2)

# spatstat object ----

min.lon <- min(V$lon)
min.lat <- min(V$lat)
V$lon <- V$lon - min.lon
V$lat <- V$lat - min.lat

# create linnet object and rescale with units meters

P <- ppp(x = V$lon, y = V$lat, 
         window = owin(xrange = c(min(V$lon), max(V$lon)), yrange = c(min(V$lat), max(V$lat))))
L <- linnet(vertices = P, edges = as.matrix(E[, 2:3]))
s <- L$dpath[191, 193]/2.2
L <- spatstat::rescale(L, s = s, unitname = "Kilometers")

data.highways <- mutate(data.highways, Longitude = (Longitude - min.lon)/s,
                             Latitude = (Latitude - min.lat)/s)
L.ppp <- ppp(x = data.highways$Longitude, y = data.highways$Latitude, L$window)
L.psp <- as.psp(L)
projection <- project2segment(L.ppp, L.psp)


# choose accidents which are less than 20 meters away from L
seg <- projection$mapXY[which(projection$d < 20)]
tp <- projection$tp[which(projection$d < 20)]

L.lpp <- as.lpp(seg = seg, tp = tp, L = L)
L <- as.linnet(L.lpp)

W <- owin(xrange = c(9, 18), yrange = c(0, 8), unitname = "Kilometers")
L <- L[W]
L.lpp <- L.lpp[W]
plot(L, box = TRUE)
plot(L.lpp, box = TRUE)

delta <- 0.05
h <- 0.025
r <- 2
L <- augment.linnet(L, delta, h, r)
L.lpp <- lpp(as.ppp(L.lpp), L)

intens <- intensity.pspline.lpp(L.lpp)
sigma <- bw.lppl(L.lpp, srange = c(1, 10))
