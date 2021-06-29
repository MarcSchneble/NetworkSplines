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
library(mgcv)

# load functions
source(file = "Functions.R")
#register_google("")

# get data on car crashes
data <- read.csv(file = "Data/CarCrashesMontgomery_Incidents.csv") %>% as_tibble() %>%
  mutate(lat = as.numeric(substring(sub("\\,.*", "", Location), first = 2)),
         lon = as.numeric(gsub("[)]" ,"" , sub("^\\S+", "", Location))),
         date = as.POSIXct(substring(Crash.Date.Time, 1, 19), format = "%m/%d/%Y %H:%M:%OS", tz = "America/New_York"),
         date = date + hours(12)*(substring(Crash.Date.Time, 21, 22) == "PM")) %>%
  filter(year(date) != 2020,
         Route.Type %in% c("Maryland (State)", "Interstate (State)", "US (State)"),
         hour(date) >= 6, hour(date) <= 21) 

# construct geometric network ----

V <- read_excel("Data/VerticesMontgomery.xlsx") %>% filter(!is.na(lat), is.na(remove)) %>%
  mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>% arrange(Id)
E <- read_excel("Data/EdgesMontgomery.xlsx") %>% filter(!is.na(from))

E$from.lat <- V$lat[match(E$from, V$Id)]
E$from.lon <- V$lon[match(E$from, V$Id)]
E$to.lat <- V$lat[match(E$to, V$Id)]
E$to.lon <- V$lon[match(E$to, V$Id)]

min.lon <- min(V$lon)
min.lat <- min(V$lat)
V$lon <- V$lon - min.lon
V$lat <- V$lat - min.lat

P <- ppp(x = V$lon, y = V$lat, 
         window = owin(xrange = c(min(V$lon), max(V$lon)), yrange = c(min(V$lat), max(V$lat))))
L <- linnet(vertices = P, edges = as.matrix(E[, 2:3]))
s <- L$dpath[191, 193]/2.2
L <- spatstat.geom::rescale(L, s = s, unitname = "Kilometers")

data <- data %>%
  mutate(lon.net = (lon - min.lon)/s, lat.net = (lat - min.lat)/s)

L.ppp <- ppp(x = data$lon.net, y = data$lat.net, L$window, marks = data$Report.Number)
L.psp <- as.psp(L)
projection <- project2segment(L.ppp, L.psp)

retain <- which(projection$d < 0.1)
seg <- projection$mapXY[retain]
tp <- projection$tp[retain]
marks <- projection$Xproj$marks[retain]
L.lpp <- as.lpp(seg = seg, tp = tp, L = L, marks = marks)

W <- owin(xrange = c(15, 33), yrange = c(0, 12), unitname = "Kilometers")
L <- L[W, snip = FALSE]
L$ind.edges <-  which((E$from.lat - min.lat)/s <= 12 & (E$from.lat - min.lat)/s >= 0 &
                        (E$to.lat - min.lat)/s <= 12 & (E$to.lat - min.lat)/s >= 0 &
                        (E$from.lon - min.lon)/s <= 33 & (E$from.lon - min.lon)/s >= 15 &
                        (E$to.lon - min.lon)/s <= 33 & (E$to.lon - min.lon)/s >= 15)
X <- as.ppp(L.lpp[W, snip = FALSE])
data <- data %>% mutate(on = factor(1*(Report.Number %in% X$marks)))
plot(L, box = TRUE)

delta <- 0.05
h <- 0.025
r <- 2
L <- augment.linnet(L, delta, h, r)
L.lpp <- lpp(X, L)
plot(unmark(L.lpp))

L.lpp$data$hour <- hour(data %>% filter(on == 1) %>% pull(date)) + floor(minute(data %>% filter(on == 1) %>% pull(date))/15)/4
type <- direction <- rep(NA, sum(L$N.m))
ind <- 1
for (m in 1:L$M) {
  type[ind:(ind + L$N.m[m] - 1)] <- E$type[L$ind.edges[m]]
  direction[ind:(ind + L$N.m[m] - 1)] <- E$direction[L$ind.edges[m]]
  ind <- ind + L$N.m[m]
}
L.lpp$domain$routetype <- factor(type, levels = c("state", "interstate", "US"))
L.lpp$domain$direction <- factor(direction, levels = c("SN", "EW", "SENW", "SWNE"))


# intensity fitting ----
intens <- intensity.pspline.lpp(L.lpp) 
intens.covariates <- intensity.pspline.lpp(L.lpp, lins = c("routetype", "direction", "dist2V"), smooths = "hour") 
plot(intens, log = TRUE)
plot(intens.covariates, log = TRUE)
plot(intens, style = "width")

# plot intensity without and with covariates
pdf(file = "Plots/Montgomery_intensity.pdf", width = 10, height = 6)
par(mar=c(0, 0, 0, 1), cex = 1.6)
plot(intens/4, log = TRUE, main = "")
dev.off()

pdf(file = "Plots/Montgomery_intensity_covariates.pdf", width = 10, height = 6)
par(mar=c(0, 0, 0, 1), cex = 1.6)
plot(intens.covariates, log = TRUE, main = "")
dev.off()

# plot smooth effects
g <- ggplot(intens.covariates$effects$smooth$hour) + 
  geom_ribbon(aes(x = x, ymin = exp(lwr), ymax = exp(upr)), fill = "grey50") + 
  geom_line(aes(x = x, y = exp(y)), color = "red") + 
  theme_bw() + 
  labs(x = "Hour of the day", y = "exp(s(t))") + 
  scale_x_continuous(breaks = seq(6, 22, 2)) + 
  scale_y_continuous(breaks = seq(0.4, 1.6, 0.2))
pdf(file = "Plots/Montgomery_smooth_t.pdf", width = 6, height = 4)
print(g)
dev.off()

sigma <- bw.lppl(L.lpp)
intens.kernel <- density.lpp(L.lpp, sigma = as.numeric(sigma), dimyx = c(256, 256))
plot(intens.kernel, log = TRUE)

sigma <- bw.scott(L.lpp)
intens.kernel2d <- density.lpp(L.lpp, sigma = sigma, distance = "euclidean", dimyx = c(256, 256))
plot(intens.kernel2d, log = TRUE)

# plot network on map ----

left <- -77.22
bottom <- 38.948 
right <- -76.96
top <- 39.122

g <- ggmap(get_map(location = c(left = left, bottom = bottom, right = right, top = top), maptype = "roadmap", scale = 2)) +
  geom_segment(data = E, aes(x = from.lon, y = from.lat, xend = to.lon, yend = to.lat), size = 1.2) +
  geom_point(data = data, aes(x = lon, y = lat, color = on)) + 
  labs(x = "Longitude", y = "Latitude")

g <- ggmap(get_map(location = c(left = left, bottom = bottom, right = right, top = top), maptype = "roadmap", scale = 2)) +
  geom_segment(data = E %>% slice(L$ind.edges), aes(x = from.lon, y = from.lat, xend = to.lon, yend = to.lat), size = 1.5) +
  geom_point(data = data %>% filter(on == 1), aes(x = lon, y = lat), col = "red", alpha = 0.1) + 
  labs(x = "Longitude", y = "Latitude")

pdf(file = "Plots/MontgomeryNetwork.pdf", height = 8, width = 10)
print(g)
dev.off()
