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
  mutate(Latitude = as.numeric(substring(sub("\\,.*", "", Location), first = 2)),
         Longitude = as.numeric(gsub("[)]" ,"" , sub("^\\S+", "", Location))))

data.maryland <- filter(data, Route.Type == "Maryland (State)")
data.US.interstate <- filter(data, 
                             Route.Type %in% c("Maryland (State)", "Interstate (State)", "US (State)"))

data2 <- filter(data, Route.Type %in% "Interstate (State)")

register_google("AIzaSyDX0CVJsIDBJF8NVFPAH84oLWPfvPa335Y")
g <- ggmap(get_map(location = c(-77.08, 39.032), zoom = 12, maptype = "roadmap"), 
           xlab = "Longitude", ylab = "Latitude") +
  geom_point(data = data2, aes(x = Longitude, y = Latitude))


