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

load_GBIF_data <- function(download_key) {
  
  destfile="./data/GBIF_data.zip"
  
  # If the data has not been downloaded already, ask for credentioals and do so now
  if (!file.exists(destfile)) {
    options(gbif_user=rstudioapi::showPrompt(title = "GBIF username", message = "GBIF username", default = ""))
    options(gbif_email=rstudioapi::showPrompt(title = "GBIF e-mail address", message = "GBIF e-mail address", default = ""))
    options(gbif_pwd=rstudioapi::askForPassword("GBIF password"))
    
    download_GBIF_API(download_key=download_key,destfile_name=destfile,n_try=20,Sys.sleep_duration=30)
  }

  gbif_data <- read.table(unzip(destfile,files="occurrence.txt"),header=T,sep="\t")[c('gbifID', 'license','year','month','day','decimalLatitude', 'decimalLongitude','coordinateUncertaintyInMeters', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'specificEpithet', 'species', 'taxonRank')]
  sp::coordinates(gbif_data) <- ~decimalLongitude + decimalLatitude
  sp::proj4string(gbif_data) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

  return(gbif_data)
}

get_grid <- function(cellsize_km, grid_extent) {
  grid_extent <- extent(bbox(SpatialPoints(matrix(grid_extent, ncol = 2, byrow = TRUE), proj4string = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))))
  ncols = ceiling((grid_extent@xmax - grid_extent@xmin) / (cellsize_km * 1000))
  nrows = ceiling((grid_extent@ymax - grid_extent@ymin) / (cellsize_km * 1000))
  xmin <- grid_extent@xmin
  ymin <- grid_extent@ymin
  xmax <- xmin + 1000 * cellsize_km * nrows
  ymax <- ymin + 1000 * cellsize_km * ncols
  
  return(raster(xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, ncols=ncols, nrows=nrows, crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))

}
