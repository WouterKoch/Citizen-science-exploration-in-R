get_GBIF_credentials <- function() {
    options(gbif_user=rstudioapi::showPrompt(title = "GBIF username", message = "GBIF username", default = ""))
    options(gbif_email=rstudioapi::showPrompt(title = "GBIF e-mail address", message = "GBIF e-mail address", default = ""))
    options(gbif_pwd=rstudioapi::askForPassword("GBIF password"))
}

download_GBIF_API <- function(download_key,n_try,Sys.sleep_duration,destfile_name){
  start_time <- Sys.time()
  n_try_count <- 1
  
  download_url <- paste("http://api.gbif.org/v1/occurrence/download/request/",
                        download_key[1],sep="")
  try_download <- try(download.file(url=download_url,destfile=destfile_name,
                                    quiet=TRUE, mode="wb"),silent = TRUE)
  
  while (inherits(try_download, "try-error") & n_try_count < n_try) {   
    Sys.sleep(Sys.sleep_duration)
    n_try_count <- n_try_count+1
    try_download <- try(download.file(url=download_url,destfile=destfile_name,
                                      quiet=TRUE),silent = TRUE)
    print(paste("trying... Download link not ready. Time elapsed (min):",
                round(as.numeric(paste(difftime(Sys.time(),start_time, units = "mins"))),2)))
  }
  
  if (!file.exists(destfile_name)) {
    print(paste("Download failed after", n_try_count, "of", n_try, "tries (taking",  round(as.numeric(paste(difftime(Sys.time(),start_time, units = "mins"))),2), "minutes)"))
  }
  else {
    print(paste("Download successful after", n_try_count, "of", n_try, "tries (taking",  round(as.numeric(paste(difftime(Sys.time(),start_time, units = "mins"))),2), "minutes)"))
  }
}

load_GBIF_data <- function(download_filter) {
  
  destfile = "./data/GBIF_data.zip"
  filter_file = "./data/GBIF_filter.txt"
  file_age = as.numeric((difftime(as.Date(file.info(destfile)$ctime), Sys.Date(), units = "days")))
  
  
  # If the GBIF data has been deleted, is older than 10 days or the filter has been changed or deleted, download new data
  if(!file.exists(destfile) || !file.exists(filter_file) || file_age > 10 || paste(download_filter, collapse=', ') != paste(readChar(filter_file, file.info(filter_file)$size))) {
    cat(paste(gbif_filter, collapse=', '), file=filter_file, sep="")
    
    # Get the credentials if needed
    if(!(exists("gbif_user") && exists("gbif_email") && exists("gbif_pwd"))) {
      get_GBIF_credentials()
    }
    
    download_key <- do.call(occ_download, download_filter) %>% occ_download_meta
    download_GBIF_API(download_key=download_key,destfile_name=destfile,n_try=40,Sys.sleep_duration=30)
  }

  gbif_data <- read.table(unzip(destfile,files="occurrence.txt"), header=T, sep="\t", quote="", fill=FALSE)[c('gbifID', 'recordedBy', 'license','year','month','day','decimalLatitude', 'decimalLongitude','coordinateUncertaintyInMeters', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specificEpithet', 'species', 'taxonRank')]

  gbif_data$season <- date_to_season(paste0(gbif_data$month, "-", gbif_data$day))
  gbif_data$season <- factor(gbif_data$season, levels=c("Winter", "Spring", "Summer","Fall"), ordered=TRUE)
  
  
  # gbif_data <- gbif_data %>% dplyr::mutate(recordedBy = strsplit(as.character(recordedBy), ",")) %>% tidyr::unnest(recordedBy)
  
  
  sp::coordinates(gbif_data) <- ~decimalLongitude + decimalLatitude
  sp::proj4string(gbif_data) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

  return(gbif_data)
}

load_GeoStat_data <- function(country) {
  
  population_file = "./data/GeoStat_data.zip"
  population_url = "http://ec.europa.eu/eurostat/cache/GISCO/geodatafiles/GEOSTAT-grid-POP-1K-2011-V2-0-1.zip"
  
  if (!file.exists(population_file)) {
    download.file(url=population_url,destfile=population_file, quiet=TRUE, mode="wb")
  }
  
  geostat_data <- read.table(unzip(population_file,files="Version 2_0_1/GEOSTAT_grid_POP_1K_2011_V2_0_1.csv"),header=T, sep=",")[c('TOT_P', 'GRD_ID','CNTR_CODE')]
  geostat_data <- subset(geostat_data, CNTR_CODE==country)

  geostat_data$y <- strtoi(substr(geostat_data$GRD_ID, 5, 8)) * 1000
  geostat_data$x <- strtoi(substr(geostat_data$GRD_ID, 10, 13)) * 1000
  
  sp::coordinates(geostat_data) <- ~x + y 
  sp::proj4string(geostat_data) <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
  return(geostat_data) 
}
  
get_grid <- function(cellsize_km, grid_extent) {
  grid_extent <- extent(bbox(SpatialPoints(matrix(grid_extent, ncol = 2, byrow = TRUE), proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))))
  ncols = ceiling((grid_extent@xmax - grid_extent@xmin) / (cellsize_km * 1000))
  nrows = ceiling((grid_extent@ymax - grid_extent@ymin) / (cellsize_km * 1000))
  xmin <- grid_extent@xmin
  ymin <- grid_extent@ymin
  xmax <- xmin + 1000 * cellsize_km * ncols
  ymax <- ymin + 1000 * cellsize_km * nrows
  
  return(raster(xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, ncols=ncols, nrows=nrows, crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

}

date_to_season <- function(date_string) {
  d <- as.Date(paste0("2012-",date_string), format = "%Y-%m-%d")

    return(
    ifelse(d >= winter | d < spring, "Winter",
      ifelse(d >= spring & d < summer, "Spring",
        ifelse(d >= summer & d < fall, "Summer",
          ifelse(d >= fall & d < winter, "Fall", "Invalid")
        )
      )
    )
  )
}
