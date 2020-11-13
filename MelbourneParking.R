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
library(mgcv)
library(ggplot2)

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
delta <- 25
h <- 5
r <- 2
L <- augment.linnet(L, delta, h, r)

# read parking data
data.parking <- readRDS("Data/data_2019_clean.rds")
data.lots <- readRDS("Data/Parking_Lots.rds") %>% filter(sensor == 1) %>%
  select(lon, lat, StreetMarker) %>%
  mutate(lon = (lon - min.lon)/s, lat = (lat - min.lat)/s)
data.parking <- left_join(data.parking, data.lots, by = "StreetMarker") %>% 
  filter(h.start >= 8, h.start < 20, State == 1) %>%
  mutate(t = 0.25*floor(m.start/15),
         weekday2 = factor(weekday, levels = c("weekday", "Sat", "Sun"))) %>%
  replace_na(replace = list(weekday2 = "weekday"))

# remove parking lots with too few events
freq <- as.data.frame(table(data.parking$StreetMarker))
lots.retain <- filter(freq, Freq >= 1000) %>% pull(Var1) 
data.parking <- filter(data.parking, StreetMarker %in% lots.retain)
data.parking$StreetMarker <- factor(data.parking$StreetMarker, levels = unique(data.parking$StreetMarker))
data.lots <- filter(data.lots, StreetMarker %in% data.parking$StreetMarker)


L.ppp <- ppp(x = data.lots$lon , y = data.lots$lat, L$window, marks = data.lots$StreetMarker)
L.psp <- as.psp(L)
projection <- project2segment(L.ppp, L.psp)

# choose parking lots which are less than 20 meters away from L
seg <- projection$mapXY[which(projection$d < 20)]
tp <- projection$tp[which(projection$d < 20)]
marks <- projection$Xproj$marks[which(projection$d < 20)]

# create L.lpp object with nearby parking lots
L.lpp <- as.lpp(seg = seg, tp = tp, L = L, marks = marks)
L <- as.linnet(L.lpp)

# intensity estimate of parking lots
intens.lots <- intensity.pspline.lpp(L.lpp)
intens.lots.covariates <- intensity.pspline.lpp(L.lpp, lins = "dist2Vdiscrete")
#intens.lots.voronoi <- densityVoronoi(L.lpp, f = 0.5, nrep = 100, dimyx = c(256, 256))

# only where intensity is at least 0.1
intens.lots.adj <- intens.lots
intens.lots.covariates.adj <- intens.lots.covariates
intens.lots.adj$v[which(intens.lots$v < 0.05, arr.ind = TRUE)] = NA
intens.lots.covariates.adj$v[which(intens.lots.covariates$v < 0.05, arr.ind = TRUE)] = NA

# plot parking lots on the geometric network
pdf(file = "Plots/MelbourneLots.pdf", width = 10, height = 8)
par(mar=c(0, 0, 0, 0), cex = 1.6)
plot(L.lpp, use.marks = FALSE, pch = 16, cols = "red", lwd = 3, legend = FALSE,
     main = "")
dev.off()

# plot intensity of parking lots
min.intens <- min(intens.lots$v, na.rm = TRUE)
max.intens <- max(intens.lots$v, na.rm = TRUE)
pdf(file = "Plots/MelbourneIntensityLots.pdf", width = 10, height = 8)
par(mar=c(0, 0, 0, 1), cex = 1.6)
plot(intens.lots.adj, log = TRUE,
     main = "")
dev.off()

pdf(file = "Plots/MelbourneIntensityLotsCovariates.pdf", width = 10, height = 8)
par(mar=c(0, 0, 0, 1), cex = 1.6)
plot(intens.lots.covariates.adj, log = TRUE,
     main = "")
dev.off()


# parking data ----

# only consider parking lots which are located on the map
data.CBD <- filter(data.parking, data.parking$StreetMarker %in% L.lpp$data$marks)
covariates <- select(data.CBD, t, weekday2)

# observed processes
L.lpp.parking <- as.lpp(x = data.CBD$lon, y = data.CBD$lat, L = L)
L.lpp.parking$data$t <- covariates$t
L.lpp.parking$data$weekday <- covariates$weekday2

# fitting
intens.parking <- intensity.pspline.lpp(L.lpp.parking, lins = "weekday", smooths = "t", eps.rho = 1e-3)

# plot smooth effects
g <- ggplot(intens.parking$effects$smooth$t) + 
  geom_line(aes(x = x, y = y)) + 
  geom_ribbon(aes(x = x, ymin = lwr, ymax = upr), color = "red") + 
  scale_x_continuous(limits = c(8, 20), breaks = 8:20) +
  theme_bw() + 
  labs(x = "Time of the day", y = "s(t)")
pdf(file = "Plots/Melbourne_smooth_t.pdf", width = 6, height = 4)
print(g)
dev.off()

intens.parking.adj <- intens.parking
intens.parking.adj$v[which(is.na(intens.lots.adj$v) | intens.parking$v < 0.2)] <- NA
intens.parking.adj$v <- intens.parking.adj$v/365*4
plot(intens.parking.adj)
plot(intens.parking.adj, log = TRUE)

intens.fluctuation <- intens.parking.adj
intens.fluctuation$v <- intens.fluctuation$v/intens.lots.adj$v

pdf("Plots/Melbourne_fluctuation_rate.pdf", width = 10, height = 6)
par(mar=c(0, 0, 0, 1), cex = 1.6)
plot(intens.fluctuation, log = TRUE, main = "")
dev.off()
