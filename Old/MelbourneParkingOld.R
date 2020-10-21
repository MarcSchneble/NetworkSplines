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

# load functions
source(file = "Functions.R")

# if FALSE load prepared RData file
load.parking.data = FALSE
edit.parking.data = FALSE


#### geometric network ####


# edges and coordinates of vertices
E = read.csv(file = "Data/EdgesWest.csv", header = TRUE, sep = ";")
V = read.csv(file = "Data/VerticesWest.csv", header = TRUE, sep = ";")

# change coordinate system
min.lon = min(V$Lon)
min.lat = min(V$Lat)
V$Lon = V$Lon - min.lon
V$Lat = V$Lat - min.lat

# create linnet object and rescale with units meters
P = ppp(x = V[, 3], y = V[, 2], 
        window = owin(xrange = c(min(V[, 3]), max(V[, 3])), yrange = c(min(V[, 2]), max(V[, 2]))))
L = linnet(vertices = P, edges = as.matrix(E[, 2:3]))
s = L$dpath[6, 26]/700
L = spatstat::rescale(L, s = s, unitname = "Meters")

# choose parameters and augment L
delta = 20
h = 1
r = 2
L = augment.linnet(L, delta, h)


#### parking lots on the geometric network ####


# read Melbourne parking data
# download: https://data.melbourne.vic.gov.au/Transport/On-street-Car-Parking-Sensor-Data-2017/u9sa-j86i

if (load.parking.data){
  Parking.Melbourne.2017 = read.csv("Data/On-street_Car_Parking_Sensor_Data_-_2017.csv", header = TRUE)
  save(Parking.Melbourne.2017, file = "RData/ParkingMelbourne2017.RData")  
} else {
  load("RData/ParkingMelbourne2017.RData")
}
Parking.Melbourne.2017 = as_tibble(Parking.Melbourne.2017)


# Melbourne parking lots in the CBD
# download: https://data.melbourne.vic.gov.au/Transport/On-street-Parking-Bays/crvt-b4kt
parking.lots = read.csv("Data/layer_0.csv", header = TRUE)

# coordinates and info whether sensor data from 2017 is available
street.markers = tibble(Marker = unique(as.vector(Parking.Melbourne.2017$StreetMarker)), Longitude = NA, Latitude = NA)
for (marker in street.markers$Marker) {
  if (is.element(marker, parking.lots$marker_id)){
    i = which(parking.lots$marker_id == marker)
    polyg = substring(gsub('.{3}$', '', parking.lots$the_geom[i]), first = 17)
    coord = strsplit(polyg, ",")[[1]]
    lon.lat = c(0, 0)
    for (j in 1:length(coord)) {
      coord[j] = trimws(gsub("\\(|\\)", "", coord[j]), which = "both")
      lon.lat = lon.lat + as.numeric(strsplit(coord[j], " +")[[1]])
    }
    lon.lat = lon.lat/length(coord)
    street.markers$Longitude[which(street.markers$Marker == marker)] = (lon.lat[1] - min.lon)/s
    street.markers$Latitude[which(street.markers$Marker == marker)] = (lon.lat[2] - min.lat)/s
  }
}

# projection of parking lots onto the geometric network
L.ppp = ppp(x = street.markers$Longitude, y = street.markers$Latitude, L$window, marks = street.markers$Marker)
L.psp = as.psp(L)
projection = project2segment(L.ppp, L.psp)

# choose parking lots which are less than 20 meters away from L
seg = projection$mapXY[which(projection$d < 20)]
tp = projection$tp[which(projection$d < 20)]
marks.western.CBD = projection$Xproj$marks[which(projection$d < 20)]

# create L.lpp object with nearby parking lots
L.lpp = as.lpp(seg = seg, tp = tp, L = L, marks = marks.western.CBD)

# intensity estimate of parking lots
intens.lots = intensity.pspline.lpp(L.lpp, r)

# only where intensity is at least 0.1
intens.lots.adj = intens.lots
intens.lots.adj$v[which(intens.lots$v < 0.1, arr.ind = TRUE)] = NA

# plot parking lots on the geometric network
pdf(file = "Plots/MelbourneLots.pdf", width = 10, height = 8)
par(mar=c(0, 1, 2, 1), cex = 2.4)
plot(L.lpp, use.marks = FALSE, pch = 16, cols = "red", lwd = 3, legend = FALSE,
     main = "On-Street Parking Lots in the \n Western CBD of Melbourne, Australia")
dev.off()

# plot intensity of parking lots
min.intens = min(intens.lots$v[-which(is.na(intens.lots$v))])
max.intens = max(intens.lots$v[-which(is.na(intens.lots$v))])
pdf(file = "Plots/MelbourneIntensityLots.pdf", width = 10, height = 8)
par(mar=c(0, 1, 2, 1), cex = 2.4)
plot(intens.lots, zlim = c(min.intens, max.intens), box = TRUE,
     main = "Intensity of On-Street Parking Lots in \n the Western CBD of Melbourne, Australia")
dev.off()


#### edit parking lot data ####


if (edit.parking.data){
  # only for the considered subset
  western.CBD = filter(Parking.Melbourne.2017, is.element(Parking.Melbourne.2017$StreetMarker, marks.western.CBD)) 
  western.CBD = as_tibble(cbind(western.CBD, street.markers[match(western.CBD$StreetMarker, street.markers$Marker), 2:3]))
  
  # change format of starting and stopping time of parking events
  kstart = strptime(as.character(western.CBD$ArrivalTime), format="%m/%d/%Y %I:%M:%S %p")
  kstop = strptime(as.character(western.CBD$DepartureTime), format="%m/%d/%Y %I:%M:%S %p")
  
  # clean up time information (not every information is needed here)
  western.CBD.clean = select(western.CBD, StreetMarker, Longitude, Latitude)
  western.CBD.clean = mutate(western.CBD.clean,
                             etype = 1*(western.CBD$Vehicle.Present == "True"),
                             enum = 0,
                             dstart = yday(strptime(as.character(western.CBD$ArrivalTime), format="%m/%d/%Y %I:%M:%S %p")),
                             dstop = yday(strptime(as.character(western.CBD$DepartureTime), format="%m/%d/%Y %I:%M:%S %p")),
                             estart.day = as.numeric(substring(kstart, 18, 19)) + 60*as.numeric(substring(kstart, 15, 16)) + 
                               3600*as.numeric(substring(kstart, 12, 13)),
                             estop.day = as.numeric(substring(kstop, 18, 19)) + 60*as.numeric(substring(kstop, 15, 16)) + 
                               3600*as.numeric(substring(kstop, 12, 13)) + 86400*(dstop == dstart+1),
                             estart = as.numeric(substring(kstart, 18, 19)) + 60*as.numeric(substring(kstart, 15, 16)) + 
                               3600*as.numeric(substring(kstart, 12, 13)) + 86400*(dstart-1),
                             estop = as.numeric(substring(kstop, 18, 19)) + 60*as.numeric(substring(kstop, 15, 16)) + 
                               3600*as.numeric(substring(kstop, 12, 13)) + 86400*(dstop-1),
                             gstart = 0,
                             gstop = western.CBD$DurationSeconds,
                             gtime = gstop,
                             gtime2 = as.double(kstop - kstart)
  )
  western.CBD.clean = arrange(western.CBD.clean, estart)
  
  # correct two consecutive events of the same type
  western.CBD.clean = mutate(western.CBD.clean, delete = 0)
  western.CBD.clean.new = filter(western.CBD.clean, delete == 1)
  
  for(marker in marks){
    subset = filter(western.CBD.clean, StreetMarker == marker)
    
    ind_delete = c(which(subset$etype[1:(nrow(subset)-1)] == subset$etype[2:nrow(subset)] & 
                           subset$estop[1:(nrow(subset)-1)] == subset$estart[2:nrow(subset)] & 
                           subset$StreetMarker[1:(nrow(subset)-1)] == subset$StreetMarker[2:nrow(subset)])+1, 1)
    subset$delete[ind_delete] = 1
    
    i = 1
    while(i < length(ind_delete)) {
      j = 1
      while (ind_delete[i] == (ind_delete[i+j]-j)) {
        j = j+1
      }
      subset$dstop[ind_delete[i]-1] = subset$dstop[ind_delete[i]-1+j]
      subset$estop.day[ind_delete[i]-1] =  subset$estop.day[ind_delete[i]-1+j]
      subset$estop[ind_delete[i]-1] =  subset$estop[ind_delete[i]-1+j]
      subset$gstop[ind_delete[i]-1] = sum(subset$gstop[(ind_delete[i]-1):(ind_delete[i]-1+j)]) 
      subset$gtime[ind_delete[i]-1] = sum(subset$gtime[(ind_delete[i]-1):(ind_delete[i]-1+j)])  
      subset$gtime2[ind_delete[i]-1] = sum(subset$gtime2[(ind_delete[i]-1):(ind_delete[i]-1+j)])  
      i = i+j
    }
    
    subset = filter(subset, delete == 0)  
    subset = select(subset, -delete)
    
    western.CBD.clean.new = rbind(western.CBD.clean.new, subset)
  }
  western.CBD.clean = western.CBD.clean.new
  
  # add weekday information
  western.CBD.clean = mutate(western.CBD.clean, weekd = (dstart-1)%%7 + 7*((dstart-1)%%7 == 0))
  
  save(western.CBD.clean, file = "RData/ParkingWesternCBD.RData")
} else {
  load("RData/ParkingWesternCBD.RData")
}


#### fit intensities ####


# extract clearing events in the morning between 8 am and 9 am / 5 pm and 6 pm
events.morning = filter(western.CBD.clean, weekd <= 5 & etype == 0 & estart.day >= 8*3600 & estart.day <= 9*3600)
events.evening = filter(western.CBD.clean, weekd <= 5 & etype == 0 & estart.day >= 17*3600 & estart.day <= 18*3600)

# observed processes
L.lpp.morning = as.lpp(x = events.morning$Longitude, y = events.morning$Latitude, L = L)
L.lpp.evening = as.lpp(x = events.evening$Longitude, y = events.evening$Latitude, L = L)

# intensity (per hour) in the morning and in the evening
intens.morning = intensity.pspline.lpp(L.lpp.morning, r)/260
intens.evening = intensity.pspline.lpp(L.lpp.evening, r)/260

# fluctuation rates 
intens.morning.normed = intens.morning
intens.morning.normed$v = intens.morning$v/intens.lots.adj$v
intens.evening.normed = intens.evening
intens.evening.normed$v = intens.evening$v/intens.lots.adj$v

# plot fluctuation rates where parking lot intensity >= 0.1
min.intens = min(min(intens.morning.normed$v[-which(is.na(intens.morning.normed$v))]), 
                 min(intens.evening.normed$v[-which(is.na(intens.evening.normed$v))]))
max.intens = max(max(intens.morning.normed$v[-which(is.na(intens.morning.normed$v))]), 
                 max(intens.evening.normed$v[-which(is.na(intens.evening.normed$v))]))

pdf(file = "Plots/MelbourneIntensityMorning.pdf", width = 10, height = 8)
par(mar=c(0, 1, 2, 1), cex = 2.4)
plot.linim(intens.morning.normed, main = "Fluctuation Rate of On-Street Parking Lots \n between 8 and 9 am", 
           zlim = c(min.intens, max.intens), box = TRUE)
dev.off()

pdf(file = "Plots/MelbourneIntensityEvening.pdf", width = 10, height = 8)
par(mar=c(0, 1, 2, 1), cex = 2.4)
plot.linim(intens.evening.normed, main = "Fluctuation Rate of On-Street Parking Lots \n between 5 and 6 pm", 
           zlim = c(min.intens, max.intens), box = TRUE) 
dev.off()
