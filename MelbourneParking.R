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
L <- spatstat.geom::rescale(L, s = s, unitname = "Meters")

# choose parameters and augment L
delta <- 10
h <- 2
r <- 2
L <- augment.linnet(L, delta, h, r)

# read parking data
data.parking <- readRDS("Data/data_2019_clean.rds")
data.lots <- readRDS("Data/Parking_Lots.rds") %>% filter(sensor == 1) %>%
  dplyr::select(lon, lat, StreetMarker) %>%
  mutate(lon = (lon - min.lon)/s, lat = (lat - min.lat)/s)
data.parking <- left_join(data.parking, data.lots, by = "StreetMarker") %>% 
  filter(h.start >= 8, h.start < 20, State == 1, month(ArrivalTime) %in% c(6, 7, 8)) %>%
  mutate(t = 0.25*floor(m.start/15),
         weekday2 = factor(weekday, levels = c("weekday", "Sat", "Sun"))) %>%
  replace_na(replace = list(weekday2 = "weekday")) 

# remove parking lots with too few events
freq <- as.data.frame(table(data.parking$StreetMarker))
lots.retain <- filter(freq, Freq >= 1*92) %>% pull(Var1) 
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
L.lpp$data$sos <- data.parking$SideOfStreetCode[match(L.lpp$data$marks, data.parking$StreetMarker)]

# intensity estimate of parking lots
intens.lots <- intensity.pspline.lpp(L.lpp, dimyx = c(150, 250))
intens.lots.covariates <- intensity.pspline.lpp(L.lpp, lins = c("dist2Vdiscrete"),
                                                dimyx = c(150, 250))

# only where intensity is at least 0.1
intens.lots.adj <- intens.lots
intens.lots.covariates.adj <- intens.lots.covariates
intens.lots.adj$v[which(intens.lots$v < 0.1, arr.ind = TRUE)] = NA
intens.lots.covariates.adj$v[which(intens.lots.covariates$v < 0.1, arr.ind = TRUE)] = NA
intens.lots.covariates.vertex <- adjust.for.vertex.distance(intens.lots.covariates, L)
intens.lots.covariates.vertex.adj <- intens.lots.covariates.vertex
intens.lots.covariates.vertex.adj$v[which(intens.lots.covariates.vertex$v < 0.1, arr.ind = TRUE)] <- NA

# plot parking lots on the geometric network
pdf(file = "Plots/MelbourneLots.pdf", width = 8, height = 5)
par(mar=c(0, 0, 0, 0), cex = 1.6)
plot(as.linnet(L.lpp), lwd = 3, 
     main = "")
points(L.lpp$data$x, L.lpp$data$y, lwd = 0.5)
dev.off()

# plot intensity of parking lots
min.intens <- min(intens.lots$v, na.rm = TRUE)
max.intens <- max(intens.lots$v, na.rm = TRUE)
pdf(file = "Plots/MelbourneIntensityLots.pdf", width = 8, height = 5)
par(mar=c(0, 0, 0, 1), cex = 1.7)
plot(intens.lots.adj, log = TRUE, box = TRUE, ribwid = 0.06, ribsep = 0.00,
     main = "")
for (e in 1:L$lines$n) {
  lines(L$lines$ends[e, c(1, 3)], L$lines$ends[e, c(2, 4)], lwd = 3, col = "grey")
}
plot(intens.lots.adj, log = TRUE,
     main = "", add = TRUE, ribwid = 0.06, ribsep = 0.00)
dev.off()

pdf(file = "Plots/MelbourneIntensityLotsCovariates.pdf", width = 8, height = 5)
par(mar=c(0, 0, 0, 1), cex = 1.7)
plot(intens.lots.covariates.adj, log = TRUE, ribwid = 0.06, ribsep = 0.00,
     main = "", box = TRUE)
for (e in 1:L$lines$n) {
  lines(L$lines$ends[e, c(1, 3)], L$lines$ends[e, c(2, 4)], lwd = 3, col = "grey")
}
plot(intens.lots.covariates.adj, log = TRUE, box = TRUE,
     main = "", add = TRUE, ribwid = 0.06, ribsep = 0.00)
dev.off()


# parking data ----

# only consider parking lots which are located on the map
data.CBD <- filter(data.parking, data.parking$StreetMarker %in% L.lpp$data$marks)
covariates <- select(data.CBD, t, weekday, SideOfStreetCode)

# observed processes
L.lpp.parking <- as.lpp(x = data.CBD$lon, y = data.CBD$lat, L = L)
L.lpp.parking$data$t <- covariates$t
L.lpp.parking$data$weekday <- covariates$weekday
L.lpp.parking$data$sos <- covariates$SideOfStreetCode

# fitting
intens.parking.covariates <- intensity.pspline.lpp(L.lpp.parking, smooths = "t",
                                                   dimyx = c(150, 250))

# plot smooth effects
g <- ggplot(intens.parking.covariates$effects$smooth$t) + 
  geom_ribbon(aes(x = x, ymin = lwr, ymax = upr), color = "grey80") + 
  geom_line(aes(x = x, y = y), color = "red") + 
  scale_x_continuous(limits = c(8, 20), breaks = 8:20) +
  theme_bw() + 
  theme(axis.title = element_text(size = 15),
        axis.text  = element_text(size = 12)) +
  labs(x = "Hour of the day", y = "s(t)")
pdf(file = "Plots/Melbourne_smooth_t.pdf", width = 6, height = 4)
print(g)
dev.off()

# fluctuation
intens.parking.adj <- intens.parking.covariates
intens.parking.adj$v[which(is.na(intens.lots.covariates.vertex.adj$v) | intens.parking.covariates$v < 1)] <- NA
intens.parking.adj$v <- intens.parking.adj$v/92*4
plot(intens.parking.adj)
plot(intens.parking.adj, log = TRUE)

intens.fluctuation <- intens.parking.adj
intens.fluctuation$v <- intens.fluctuation$v/intens.lots.covariates.vertex.adj$v

ind <- which(intens.fluctuation$v > 5, arr.ind = TRUE)
y.coord <- intens.fluctuation$yrow[ind[, 1]]
x.coord <- intens.fluctuation$xcol[ind[, 2]]

pdf("Plots/Melbourne_fluctuation_rate.pdf", width = 8, height = 5)
par(mar=c(0, 0, 0, 1), cex = 1.6)
plot(intens.fluctuation, log = TRUE, main = "", ribwid = 0.06, ribsep = 0.00, box = TRUE)
for (e in 1:L$lines$n) {
  lines(L$lines$ends[e, c(1, 3)], L$lines$ends[e, c(2, 4)], lwd = 3, col = "grey")
}
points(x.coord, y.coord, cex = 2)
plot(intens.fluctuation, log = TRUE, main = "", add = TRUE, ribwid = 0.06, ribsep = 0.00, box = TRUE)
dev.off()
