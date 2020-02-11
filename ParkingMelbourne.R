library(spatstat)
library(igraph)
library(Matrix)
library(MASS)
library(VCA)
library(psych)
library(colorRamps)
library(dplyr)
library(splines)

setwd("X:/Projekte/Network Splines/R Code/Git")
source(file = "Functions.R")

# Extract of CBD network
#V = read.csv(file = "Vertices.csv", header = TRUE, sep = ";")
#E = read.csv(file = "Edges.csv", header = TRUE, sep = ";")
#save(V, file = "VerticesCBD.RData")
#save(E, file = "EdgesCBD.RData")
load("VerticesCBD.RData")
load("EdgesCBD.RData")

# change coordinate system
min.lon = min(V$Lon)
min.lat = min(V$Lat)
V$Lon = V$Lon - min.lon
V$Lat = V$Lat - min.lat

# create linnet object and rescale with units meters
P = ppp(x = V[, 3], y = V[, 2], 
        window = owin(xrange = c(min(V[, 3]), max(V[, 3])), yrange = c(min(V[, 2]), max(V[, 2]))))
L = linnet(vertices = P, edges = as.matrix(E[, 2:3]))
s = L$dpath[1, 3]/440
L = spatstat::rescale(L, s = s, unitname = "Meters")



# read Melbourne parking data
#Data = read.csv("On-street_Car_Parking_Sensor_Data_-_2017.csv", header = TRUE)
#save(Data, file = "ParkingMelbourne2017.RData")
load("ParkingMelbourne2017.RData")

#bays = read.csv("layer_0.csv", header = TRUE)
#save(bays, file = "BaysCBD.RData")
load("BaysCBD.RData")


# extract only streets of the map extract
sub.Data = filter(Data, StreetName == "SPRING STREET" | StreetName == "EXHIBITION STREET" | StreetName == "RUSSELL STREET"
                  | StreetName == "FLINDERS STREET" | StreetName == "COLLINS STREET" | StreetName == "LITTLE COLLINS STREET"
                  | StreetName == "BOURKE STREET" | StreetName == "LITTLE BOURKE STREET" | StreetName == "LONSDALE STREET"
                  | StreetName == "LITTLE LONSDALE STREET" | StreetName == "LA TROBE STREET") 

# add coordinates to parking data
sub.Data = mutate(sub.Data, Longitude = NA, Latitude = NA)
markers = unique(as.vector(sub.Data$StreetMarker))
k = 0
for (marker in markers) {
  k = k+1
  print(k)
  if (is.element(marker, bays$marker_id)){
    i = which(bays$marker_id == marker)
    polyg = substring(gsub('.{3}$', '', bays$the_geom[i]), first = 17)
    co = strsplit(polyg, ",")[[1]]
    lon.lat = c(0, 0)
    for (j in 1:length(co)) {
      co[j] = trimws(gsub("\\(|\\)", "", co[j]), which = "both")
      lon.lat = lon.lat + as.numeric(strsplit(co[j], " +")[[1]])
    }
    lon.lat = lon.lat/length(co)
    sub.Data$Longitude[which(sub.Data$StreetMarker == marker)] = lon.lat[1] - min.lon
    sub.Data$Latitude[which(sub.Data$StreetMarker == marker)] = lon.lat[2] - min.lat
  }
}

# extract only the parking lots where coordinate info is available
sub.Data = filter(sub.Data, !is.na(Longitude))


###############################################################################################################################
###############################################################################################################################

# extract only parking lots of the map ectract


# Exhibition Street
range(sub.Data[which(sub.Data$StreetName == "EXHIBITION STREET"), ]$Latitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "EXHIBITION STREET" & sub.Data$Latitude > -37.808132 - min.lat), ]





# Russell Street
range(sub.Data[which(sub.Data$StreetName == "RUSSELL STREET"), ]$Latitude) 
sub.Data = sub.Data[-which(sub.Data$StreetName == "RUSSELL STREET" & sub.Data$Latitude > -37.808851 - min.lat), ] 


# Flinders Street
range(sub.Data[which(sub.Data$StreetName == "FLINDERS STREET"), ]$Longitude) 
sub.Data = sub.Data[-which(sub.Data$StreetName == "FLINDERS STREET" & sub.Data$Longitude < 144.969825 - min.lon), ] 


# Collins Street
range(sub.Data[which(sub.Data$StreetName == "COLLINS STREET"), ]$Longitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "COLLINS STREET" & sub.Data$Longitude < 144.968872 - min.lon), ]


# Little Collins Street
range(sub.Data[which(sub.Data$StreetName == "LITTLE COLLINS STREET"), ]$Longitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "LITTLE COLLINS STREET" & sub.Data$Longitude < 144.968441 - min.lon), ]


# Bourke Street
range(sub.Data[which(sub.Data$StreetName == "BOURKE STREET"), ]$Longitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "BOURKE STREET" & sub.Data$Longitude < 144.968013 - min.lon), ]


# Little Bourke Street
range(sub.Data[which(sub.Data$StreetName == "LITTLE BOURKE STREET"), ]$Longitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "LITTLE BOURKE STREET" & sub.Data$Longitude < 144.967567 - min.lon), ]


# Lonsdale Street
range(sub.Data[which(sub.Data$StreetName == "LONSDALE STREET"), ]$Longitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "LONSDALE STREET" & sub.Data$Longitude < 144.967085 - min.lon), ]


# Little Lonsdale Street
range(sub.Data[which(sub.Data$StreetName == "LITTLE LONSDALE STREET"), ]$Longitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "LITTLE LONSDALE STREET" & sub.Data$Longitude < 144.966616 - min.lon), ]


# La Trobe Street
range(sub.Data[which(sub.Data$StreetName == "LA TROBE STREET"), ]$Longitude)
sub.Data = sub.Data[-which(sub.Data$StreetName == "LA TROBE STREET" & sub.Data$Longitude < 144.966200 - min.lon), ]


###############################################################################################################################
###############################################################################################################################


sub.Data = as_tibble(sub.Data)
markers = as.character(unique(sub.Data$StreetMarker))
markers = markers[order(nchar(markers), markers)]

kstart = strptime(as.character(sub.Data$ArrivalTime), format="%m/%d/%Y %I:%M:%S %p")
kstop = strptime(as.character(sub.Data$DepartureTime), format="%m/%d/%Y %I:%M:%S %p")

Parking = select(sub.Data, StreetMarker, Longitude, Latitude)
Parking = mutate(Parking,
                 etype = 1*(sub.Data$Vehicle.Present == "True"),
                 enum = 0,
                 dstart = yday(strptime(as.character(sub.Data$ArrivalTime), format="%m/%d/%Y %I:%M:%S %p")),
                 dstop = yday(strptime(as.character(sub.Data$DepartureTime), format="%m/%d/%Y %I:%M:%S %p")),
                 estart_d = as.numeric(substring(kstart, 18, 19)) + 60*as.numeric(substring(kstart, 15, 16)) + 
                   3600*as.numeric(substring(kstart, 12, 13)),
                 estop_d = as.numeric(substring(kstop, 18, 19)) + 60*as.numeric(substring(kstop, 15, 16)) + 
                   3600*as.numeric(substring(kstop, 12, 13)) + 86400*(dstop == dstart+1),
                 estart = as.numeric(substring(kstart, 18, 19)) + 60*as.numeric(substring(kstart, 15, 16)) + 
                   3600*as.numeric(substring(kstart, 12, 13)) + 86400*(dstart-1),
                 estop = as.numeric(substring(kstop, 18, 19)) + 60*as.numeric(substring(kstop, 15, 16)) + 
                   3600*as.numeric(substring(kstop, 12, 13)) + 86400*(dstop-1),
                 gstart = 0,
                 gstop = sub.Data$DurationSeconds,
                 gtime = gstop,
                 gtime2 = as.double(kstop - kstart)
)
Parking = arrange(Parking, estart)


Parking = mutate(Parking, delete = 0)
Parking_new = filter(Parking, delete == 1)

for(marker in markers){
  print(marker)
  subset = filter(Parking, StreetMarker == marker)
  
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
    subset$estop_d[ind_delete[i]-1] =  subset$estop_d[ind_delete[i]-1+j]
    subset$estop[ind_delete[i]-1] =  subset$estop[ind_delete[i]-1+j]
    subset$gstop[ind_delete[i]-1] = sum(subset$gstop[(ind_delete[i]-1):(ind_delete[i]-1+j)]) 
    subset$gtime[ind_delete[i]-1] = sum(subset$gtime[(ind_delete[i]-1):(ind_delete[i]-1+j)])  
    subset$gtime2[ind_delete[i]-1] = sum(subset$gtime2[(ind_delete[i]-1):(ind_delete[i]-1+j)])  
    i = i+j
  }
  
  subset = filter(subset, delete == 0)  
  subset = select(subset, -delete)
  
  Parking_new = rbind(Parking_new, subset)
}

# Test if everything went right
ind_delete = which(Parking_new$etype[1:(nrow(Parking_new)-1)] == Parking_new$etype[2:nrow(Parking_new)] & 
                     Parking_new$estop[1:(nrow(Parking_new)-1)] == Parking_new$estart[2:nrow(Parking_new)] & 
                     Parking_new$StreetMarker[1:(nrow(Parking_new)-1)] == Parking_new$StreetMarker[2:nrow(Parking_new)])+1
ind_delete
Parking = Parking_new

# add external information
Parking = mutate(Parking, weekd = (dstart-1)%%7 + 7*((dstart-1)%%7 == 0))
Parking$StreetMarker = as.character(Parking$StreetMarker)

# extract only events that start after 8 am and end before 8 pm at the same day
select = which(Parking$dstart == Parking$dstop & Parking$estart_d > 8*3600 & Parking$estop_d < 20*3600 & Parking$gtime > 0)
Parking = Parking[select, ]


#save(Parking, file = "Parking.RData")
load("Parking.RData")
markers = unique(Parking$StreetMarker)

# average duration and number of events
average.duration = n.events = matrix(0, length(unique(Parking$StreetMarker)), 2)

for (i in 1:nrow(average.duration)) {
  # No parking event
  average.duration[i, 1] = mean(Parking$gtime[which(Parking$StreetMarker == markers[i] & 
                                                      Parking$etype == 0)])
  n.events[i, 1] = length(Parking$gtime[which(Parking$StreetMarker == markers[i] & 
                                                Parking$etype == 0)])
  # Parking event
  average.duration[i, 2] = mean(Parking$gtime[which(Parking$StreetMarker == markers[i] & 
                                                      Parking$etype == 1)])
  n.events[i, 2] = length(Parking$gtime[which(Parking$StreetMarker == markers[i] & 
                                                Parking$etype == 1)])
}
average.duration


# network
delta = 10
h = 1
r = 2
L = augment.linnet(L, delta, h)

# x and y coordinates
Parking$Longitude = Parking$Longitude/s
Parking$Latitude = Parking$Latitude/s

weekdays = filter(Parking, weekd <= 5 & etype == 0 & dstart <= 31)
weekends = filter(Parking, weekd >= 6 & etype == 0 & dstart <= 31)

L.lpp.lots = as.lpp(x = unique(Parking$Longitude), y = unique(Parking$Latitude), L = L)
L.lpp.weekdays = as.lpp(x = weekdays$Longitude, y = weekdays$Latitude, L = L)
L.lpp.weekends = as.lpp(x = weekends$Longitude, y = weekends$Latitude, L = L)



intens.weekdays = intensity.psplines.lpp(L.lpp.weekdays, r, rho = 10)/(22*12)
intens.weekends = intensity.psplines.lpp(L.lpp.weekends, r, rho = 10)/(9*12)

range.intens = range=c(min(intens.weekdays$v[-which(is.na(intens.weekdays$v))], intens.weekends[-which(is.na(intens.weekends$v))]), 
                       max(intens.weekdays$v[-which(is.na(intens.weekdays$v))], intens.weekends[-which(is.na(intens.weekends$v))]))

# plots
pdf(file = "Plots/EasternCBD.pdf", width = 10, height = 8)
par(mar=c(0, 0, 2, 0.5), cex = 2)
plot(L.lpp.lots, main = "Parking Bays with installed in-Ground Sensors \n in the Eastern CBD of Melbourne, Australia", lwd = 3, cols = "red", chars = 16, size = 1)
dev.off()

pdf(file = "Plots/IntensityWeekdays.pdf", width = 10, height = 8)
par(mar=c(0, 0, 1, 0.5), cex = 2)
plot(intens.weekdays, style = "width", zlim = range.intens, scale = 40, main = "Estimated Hourly Intensity on Weekdays")
dev.off()

pdf(file = "Plots/IntensityWeekend.pdf", width = 10, height = 8)
par(mar=c(0, 0, 1, 0.5), cex = 2)
plot(intens.weekends, style = "width", zlim = range.intens, scale = 40, main = "Estimated Hourly Intensity on the Weekend")
dev.off()
