read_gbif <- function(path) {
  result <- read.table(path, sep="\t", header=TRUE, quote = "", fill=TRUE)
  result
}

read_geostat <- function(path) {
  result <- read.table(path, sep=",", header=TRUE, quote = "", fill=TRUE)
  result
}

geostat_add_utm <- function(dataset) {
  dataset$utm_y <- strtoi(substr(dataset$X.GRD_ID, 6, 9)) * 1000
  dataset$utm_x <- strtoi(substr(dataset$X.GRD_ID, 11, 14)) * 1000
  
  xy <- data.frame(ID = 1:length(dataset[, 1]), X = dataset$utm_x, Y = dataset$utm_y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  res <- as.data.frame(res)
  
  dataset$long <- res$X
  dataset$lat <- res$Y

  return(dataset)
}

filter_by_genus <- function(dataset, genera) {
  result <- subset(dataset, genus %in% genera & taxonrank == "SPECIES")
  result <- droplevels(result)
  result
}

summarize_list <- function(list) {
  if (length(list) > 1) {
    list[length(list)-1] <- paste(list[length(list)-1], list[length(list)], sep=" and ")
    list <- list[-length(list)]
  }
  paste(list, collapse = ", ")
}

add_rounded_coordinates <- function(dataset, km=1) {
  xy <- data.frame(ID = 1:length(dataset[, 1]), X = dataset$decimallongitude, Y = dataset$decimallatitude)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
  res <- spTransform(xy, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
  res <- as.data.frame(res)
  
  dataset$utm_x <- res$X + 500 * km
  dataset$utm_y <- res$Y + 500 * km 
  dataset$utm_x_rounded <- (dataset$utm_x / 1000) %/% km
  dataset$utm_y_rounded <- (dataset$utm_y / 1000) %/% km
  
  xy <- data.frame(ID = 1:length(dataset[, 1]), X = dataset$utm_x_rounded * 1000 * km, Y = dataset$utm_y_rounded * 1000 * km)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  res <- as.data.frame(res)
  
  dataset$long_rounded <- res$X
  dataset$lat_rounded <- res$Y

  return(dataset)
}



utm_add_rasterized <- function(dataset, km=1) {
  dataset$utm_x_rounded <- ((dataset$utm_x + 500 * km) / 1000) %/% km
  dataset$utm_y_rounded <- ((dataset$utm_y + 500 * km) / 1000) %/% km
  
  dataset <- dataset %>% group_by(utm_x_rounded, utm_y_rounded) %>% summarise(pop=sum(X.TOT_P.))

  xy <- data.frame(ID = 1:length(dataset[, 1]), X = dataset$utm_x_rounded * 1000 * km, Y = dataset$utm_y_rounded * 1000 * km)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
  res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
  res <- as.data.frame(res)
  
  dataset$long_rounded <- res$X
  dataset$lat_rounded <- res$Y

  dataset$raster_cell <- paste("x", dataset$utm_x_rounded, "y", dataset$utm_y_rounded)

  return(dataset)
}

join_by_utm_rounded <- function(observations, environment) {
  observations <- observations %>% group_by(utm_x_rounded, utm_y_rounded) %>% summarise(num=n())
  return(observations %>% inner_join(environment, by = c("utm_x_rounded", "utm_y_rounded")))
}







