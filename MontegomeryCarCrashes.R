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

data <- read.csv(file = "Data/CarCrashesMontgomery.csv") %>% as_tibble() %>%
  mutate(Latitude = as.numeric(substring(sub("\\,.*", "", data$Location), first = 2)),
         Longitude = as.numeric(gsub("[)]" ,"" , sub("^\\S+", "", data$Location))))

data.maryland <- filter(data, Route.Type == "Maryland (State)")
data.US.interstate <- filter(data, 
                             Route.Type %in% c("Maryland (State)", "Interstate (State)", "US (State)", "Ramp", "Government"))

register_google("AIzaSyDX0CVJsIDBJF8NVFPAH84oLWPfvPa335Y")
g <- ggmap(get_map(location = c(-77.2, 39.14), zoom = 12, maptype = "roadmap"), 
           xlab = "Longitude", ylab = "Latitude") +
  geom_point(data = data.US.interstate, aes(x = Longitude, y = Latitude))


