#' Adds cruise behavior
#'
#' Finds cruise behavior that meet given threshold parameters
#'
#' @usage AddCruiseBehavior(df, alt_abv_grd, min_speed)
#' @param df dataframe
#' @param min_agl distance (in meters) above ground
#' @param min_speed minimum speed of bird
#'
#' @return dataframe with behavior column that has "cruise"
#' @export
#'
#' @details should run AddNestBehavior, AddRoostBehavior, and AddCruiseBehavior
#'   prior to this function
AddCruiseBehavior <- function(df = df,
                              min_agl = 200,
                              min_speed = 5) {
  df <- df
  df$behavior <- ifelse(df$agl >= min_agl & df$speed >= min_speed &
    is.na(df$behavior), "cruise", df$behavior)
  return(df)
}

#' Adds flight behavior
#'
#' Finds location data that meet given threshold parameters for flights, which
#'   include the start point, mid-flight, and final location of a flight.
#'
#' @param df dataframe
#' @param min_speed minimum speed of bird
#' @param min_step_length distance (in meters) between locations
#' @param max_step_time maximum time between locations
#'
#' @return dataframe with columns for flight_index, flight_step, and
#'  flight_length
#' @export
#'
#' @details should run AddLandscapeValues() prior to this function
AddFlightData <- function(df,
                          min_speed = 5,
                          min_step_length = 50,
                          max_step_time = 20){
  df_org <- df
  df <- df
  df$flight <- ifelse(df$speed > min_speed, TRUE, FALSE)
  df$movement <- ifelse(df$step_length > min_step_length, TRUE, FALSE)
  df$within_window <- ifelse(df$step_time < max_step_time, TRUE , FALSE)
  moved <- ifelse(df$step_length > min_step_length, TRUE, FALSE)
  df$moved <- c(NA, moved[-length(moved)])
  df$mid_flight <- ifelse(df$moved == TRUE & df$movement == TRUE &
    df$flight == TRUE, TRUE, FALSE)
  df$index <- seq.int(nrow(df))
  before <- which(df$mid_flight == TRUE) - 1
  mid <- which(df$mid_flight == TRUE)
  after <- which(df$mid_flight == TRUE) + 1
  complete <- sort.int(unique(c(before, mid, after)))
  flights <- df[complete,c("id", "datetime", "index", "mid_flight")]
  index1 <- c(which(!diff(flights$index)==1), nrow(flights))
  index2 <- c(1, index1 + 1)  # index of first records of sequential groups
  index3 <- diff(index2)  # length of each seq group
  index4 <- rep(seq(1,length(index3), by=1), times=index3)
  flights$flight_index <- index4
  flights$flight_step <- sequence(rle(flights$flight_index)$lengths)
  flights$flight_index <- factor(index4)
  flights <- plyr::ddply(flights, .(flight_index, id), transform,
    flight_length=length(flight_step))
  df2 <- suppressMessages(plyr::join(df_org, flights))
  df2$index <- NULL
  return(df2)
}

#' Adds home and conspecific edge distances to 'baea' dataframe
#'
#' Creates Raster Layers of homerange kernels based on home centroids
#'
#' @param baea dataframe of baea locations
#' @param nest_set dataframe of nests
#' @param base base Raster that sets the projection, extent, and dimensions of
#'   the study area
#' @param output_dir directory for output files (distance, homerange)
#' @param max_r maximum radius to calculate the homerange raster from each
#'   df_home centroid
#' @param write_home_dist logical, write home nest distance raster to file.
#'   Default is FALSE.
#' @param write_con_dist logical, write conspecific nest distance raster to
#'   file. Default is FALSE.
#' @param write_con_dist_nest logical, write conspecific nest distance with all
#'   values centered on zero at the nest as raster to file. Default is FALSE.
#' @param write_edge_dist logical, write territorial edge distance raster to
#'   file. Default is FALSE.
#' @param write_terr_edge logical, write territorial edge raster to file.
#'   Default is FALSE.
#'
#' @return baea
#'
#' @importFrom magrittr "%>%"
#' @export
#'
#' @details Creates folders for every year in output directory
AddHomeConEdgeDistanceBAEA <- function(baea,
                                       nest_set,
                                       base,
                                       output_dir = getwd(),
                                       max_r,
                                       write_home_dist = FALSE,
                                       write_con_dist = FALSE,
                                       write_con_dist_nest = FALSE,
                                       write_edge_dist = FALSE,
                                       write_terr_edge = FALSE){
  id <- "nest_site"
  name <- "name"
#  if (output_dir != getwd()) unlink(output_dir, recursive=TRUE)
  if (!dir.exists(output_dir)) dir.create(output_dir)
  if (write_home_dist | write_edge_dist | write_edge_dist | write_terr_edge){
    for (k in sort(unique(baea$year))) dir.create(file.path(output_dir, k),
      showWarnings = FALSE)
  }
  distance_df <- data.frame(year = numeric(), id=character(),
    home_dist_min=numeric(), home_dist_max=numeric(), con_dist_min=numeric(),
    con_dist_max=numeric(), edge_dist_min=numeric(), edge_dist_max=numeric())
  cellsize <- res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  for (i in 1:nrow(nest_set)){
    nest_set[i, "x"] <- CenterXYInCell(nest_set[i,"long_utm"],
      nest_set[i,"lat_utm"], xmin, ymin, cellsize)[1]
    nest_set[i, "y"] <- CenterXYInCell(nest_set[i,"long_utm"],
      nest_set[i,"lat_utm"], xmin, ymin, cellsize)[2]
  }
  for (i in sort(unique(baea$year))){
    baea_year <- baea %>% dplyr::filter(year == i)
    active_year <- paste0("active_", i)
    col <- which(colnames(nest_set) == active_year)
    df_all <- nest_set[which(nest_set[,col] == TRUE), ]
    year_nest <- unique(baea_year$nest_site)
    df_home <- df_all %>% dplyr::filter(nest_site %in% year_nest)
    df_all_sp <- SpatialPointsDataFrame(df_all[,c("x","y")], df_all,
      proj4string=crs(base))
    homerange_ids <- df_home[,id]
    total <- length(homerange_ids)
    for (j in 1:nrow(df_home)) {
      hr_nest_site <- df_home[j, id]
      id_year_xy <- baea %>%
        dplyr::filter(nest_site == hr_nest_site) %>%
        dplyr::filter(year == i) %>%
        dplyr::mutate(x = long_utm) %>%
        dplyr::mutate (y = lat_utm) %>%
        dplyr::select(x,y)
      sv <- baea$nest_site == hr_nest_site & baea$year == i
      sv <- ifelse(is.na(sv), FALSE, sv)
      home <- df_home[j,]
      ifelse(is.null(name), home_name <- j, home_name <- home[,name])
      writeLines(noquote(paste0("Calculating nest and edge distances for ",
        i, ": ", home_name, " (", j, " of ", total, ").")))
      home_sp <- sp::SpatialPointsDataFrame(home[,c("x","y")], home,
        proj4string=crs(base))
      xy <- c(home_sp@coords[,"x"], home_sp@coords[,"y"])
      cell_extent <- raster::extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2),
        xy[2]-(cellsize/2), xy[2]+(cellsize/2))
      cell <- setValues(raster(cell_extent,crs=projection(base),res=cellsize),j)
      home_ext <- raster::extend(cell, c(max_r_cells, max_r_cells), value=NA)
      home_dist <- raster::distanceFromPoints(home_ext, home[,c("x","y")])
      plot(home_dist)
      if (write_home_dist == TRUE) {
        filename <- file.path(output_dir, i, paste0("HomeDist_", home_name,
          ".tif"))
        raster::writeRaster(home_dist, filename=filename, format="GTiff",
          overwrite=TRUE)
        writeLines(noquote(paste("Writing:", filename)))
      }
      base_crop <- raster::raster(raster::extent(home_ext), resolution=30,
        crs=raster::crs(home_ext))
      rm(home_ext)
      base_crop <- raster::setValues(base_crop, runif(ncell(base_crop), 1, 10))
      con_dist <- raster::distanceFromPoints(base_crop,
        df_all_sp[which(df_all_sp$nest_area != home$nest_area),])
      # raster::plot(con_dist)
      if (write_con_dist == TRUE) {
        filename <- file.path(output_dir, i, paste0("ConDist_", home_name,
          ".tif"))
        writeLines(noquote(paste0("Writing: ", filename)))
        writeRaster(con_dist, filename=filename, format="GTiff",
          overwrite=TRUE)
      }
      fun <- function(x){ raster::extract(con_dist, home_sp) - x  }
      con_dist_nest <- calc(con_dist, fun)
      if (write_con_dist_nest == TRUE) {
        filename <- file.path(output_dir, i, paste0("ConDistNest_", home_name,
            ".tif"))
        writeLines(noquote(paste0("Writing: ", filename)))
        writeRaster(con_dist_nest, filename=filename, format="GTiff",
          overwrite=TRUE)
      }


      baea[sv, "con_dist"] <- raster::extract(con_dist, id_year_xy)
      baea[sv, "con_dist_min"] <- cellStats(con_dist, min)
      baea[sv, "con_dist_nest"] <- raster::extract(con_dist, home_sp)
      k <- nrow(distance_df) + 1
      distance_df[k, "con_dist_min"] <- cellStats(con_dist, min)
      distance_df[k, "con_dist_max"] <- cellStats(con_dist, max)
      distance_df[k, "con_dist_max"] <- raster::extract(con_dist, home_sp)
      rm(con_dist)
      global_dist_crop <- distanceFromPoints(base_crop, df_all_sp)
      rm(base_crop)
      cent_dist <- overlay(home_dist, global_dist_crop, fun=function (x,y)
        {ifelse(x != y, NA, x)})
      baea[sv, "home_dist"] <- raster::extract(home_dist, id_year_xy)
      distance_df[k, "home_dist_min"] <- cellStats(home_dist, min)
      distance_df[k, "home_dist_max"] <- cellStats(home_dist, max)
      cent_bounds <- boundaries(cent_dist)
      rm(home_dist)
      if (write_terr_edge == TRUE) {
        terr_edge <- cent_bounds
        terr_edge[terr_edge == 0] <- NA
        terr_edge_poly <- rasterToPolygons(terr_edge, n=4, # fun=function(x){x==1},
          digits = 8, dissolve=FALSE)
        file_dir <- file.path(output_dir, i, "TerrEdge_Shapefiles")
        if(!dir.exists(file_dir)) dir.create(file_dir)
        writeLines(noquote(paste0("Writing: ", file.path(file_dir, paste0(
          "TerrEdge_", home_name)))))
        rgdal::writeOGR(terr_edge_poly, dsn = file_dir, layer = paste0(
          "TerrEdge_", home_name), driver = "ESRI Shapefile", overwrite_layer =
          TRUE)
        rm(terr_edge_poly)
      }
      cent_bounds <- subs(cent_bounds, data.frame(from=c(0,1), to=c(NA,1)))
      edge_dist_abs <- distance(cent_bounds)
      rm(cent_bounds) # new
      edge_dist <- overlay(cent_dist, edge_dist_abs, fun=function(x,y)
        {ifelse(!is.na(x), y*-1, y*1)})
      rm(cent_dist, edge_dist_abs)
      if (write_edge_dist == TRUE) {
        filename <- file.path(output_dir, i, paste0("EdgeDist_", home_name,
          ".tif"))
        writeLines(noquote(paste0("Writing: ", filename)))
        writeRaster(edge_dist, filename=filename, format="GTiff",
          overwrite=TRUE)
        edge_dist_shift <- calc(edge_dist, function(x) x - cellStats(edge_dist,
          min))
        filename <- file.path(output_dir, i, paste0("EdgeDistShift_", home_name,
          ".tif"))
        writeLines(noquote(paste0("Writing: ", filename)))
        writeRaster(edge_dist_shift, filename=filename, format="GTiff",
          overwrite=TRUE)
        rm(edge_dist_shift)
      }
      baea[sv, "edge_dist"] <- raster::extract(edge_dist, id_year_xy)
      baea[sv, "edge_dist_min"] <- cellStats(edge_dist, min)
      rm(sv, id_year_xy)
      distance_df[k, "year"] <- i
      distance_df[k, "id"] <- home_name
      distance_df[k, "edge_dist_min"] <- cellStats(edge_dist, min)
      distance_df[k, "edge_dist_max"] <- cellStats(edge_dist, max)
      rm(edge_dist)
    }
  }
  filepath <- file.path(output_dir, "Raster_Distance_Metrics.csv")
  writeLines(noquote(paste0("Writing: ", filepath)))
  write.csv(distance_df, filepath)
  baea$con_dist_shift <- baea$con_dist - baea$con_dist_min
  baea$edge_dist_shift <- baea$edge_dist - baea$edge_dist_min
  return(baea)
}

#' Creates and exports homerange kernels
#'
#' Creates RasterLayers of homerange kernels based on home centroids
#'
#' @usage AddHomeRangeKernelsBAEA(baea, nest_set, set_name, base, output_dir,
#'   max_r, home_inflection, home_scale, avoid_inflection, avoid_scale,
#'   min_prob, write_distance, write_homerange)
#'
#' @param baea dataframe of all homerange centroids
#' @param nest_set dataframe of nests
#' @param set_name character, name of nest set
#' @param base base Raster that sets the projection, extent, and dimensions of
#'   the study area
#' @param output_dir directory for output files (distance, homerange)
#' @param max_r maximum radius to calculate the homerange raster from each
#'   df_home centroid
#' @param home_inflection inflection point of the Logistic function that governs
#'   the home kernel
#' @param home_scale scale parameter of the Logistic function that governs the
#'   scale parameter
#' @param avoid_inflection inflection point of the Logistic function that
#'   governs the conspecific avoidance kernel
#' @param avoid_scale scale parameter of the Logistic function that governs
#'   the conspecific avoidance kernel
#' @param min_prob numeric, minimum probablity threshold; all probabilities
#'   below this value will be rounded to zero. Default is 0.
#' @param write_distance logical, write distance raster to file. Default is
#'   FALSE.
#' @param write_homerange logical, write home range raster to file. Default is
#'   FALSE.
#'
#' @return A list containing homerange kernel Rasters for all the df_home
#'   centroids
#'
#' @importFrom magrittr "%>%"
#' @export
#'
#' @details Create a folder for every year in the output directory
#'
AddHomeRangeKernelsBAEA <- function(baea,
                                    nest_set,
                                    set_name,
                                    base,
                                    output_dir = getwd(),
                                    max_r,
                                    home_inflection,
                                    home_scale,
                                    avoid_inflection,
                                    avoid_scale,
                                    min_prob = 0,
                                    write_distance = TRUE,
                                    write_homerange = TRUE){
  id <- "nest_site"
  name <- "name"
#  unlink(output_dir, recursive=TRUE)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
    for (k in sort(unique(baea$year))) dir.create(file.path(output_dir, k))
  }
  for (i in sort(unique(baea$year))){
    baea_year <- baea %>% filter(year == i)
    active_year <- paste0("active_", i)
    col <- which(colnames(nest_set) == active_year)
    df_all <- nest_set[which(nest_set[,col] == TRUE), ] %>%
      mutate(x = long_utm) %>%
      mutate(y = lat_utm)
    year_nest <- unique(baea_year$nest_site)
    df_home <- df_all %>% dplyr::filter(nest_site %in% year_nest)
    homeranges_year <- CreateHomeRangeKernelsBAEA(df_all, df_home, base, max_r,
      home_inflection, home_scale, avoid_inflection, avoid_scale, min_prob,
      id, name, output_dir=file.path(output_dir, i), write_distance,
      write_homerange)
    for (j in 1:length(homeranges_year)){
      homerange_year <- homeranges_year[[j]]
      hr_nest_site <- gsub("X", "", names(homerange_year))
      id_year_xy <- baea %>%
        filter(nest_site == hr_nest_site) %>%
        filter(year == i) %>%
        mutate(x = long_utm) %>%
        mutate (y = lat_utm) %>%
        select(x,y)
      baea_name <- unique(baea %>%
        filter(nest_site == hr_nest_site) %>%
        select(id))[1,1]
      ExportKMLRasterOverlay(raster = homerange_year,
        color_pal = rev(jet2.col(255)), alpha=0.75, zip=TRUE,
        method="ngb", overwrite=TRUE, maxpixels = 300000, blur=10,
        colNA="transparent", outfile=paste0("HomeRange_", baea_name),
        output_dir=file.path(output_dir, i))
      id_year_xy_extract <- extract(homerange_year, id_year_xy)
      sv <- baea$nest_site == hr_nest_site & baea$year == i
      sv <- ifelse(is.na(sv), FALSE, sv)
      baea[sv, set_name] <- id_year_xy_extract
    }
  }
  return(baea)
}

#' Adds nest behavior to baea dataframe
#'
#' Finds nest attendance behavior that meet given threshold parameters
#'
#' @param df dataframe
#' @param distance_to_nest distance (in meters) from location to nest to assign
#'   a location as "nest" behavior
#'
#' @return dataframe with behavior column that has "nest"
#' @export
#'
#' @details need to run AddNestData() prior to running this
#'
AddNestBehavior <- function(df = df,
                            distance_to_nest = 50) {
  df<-df
  df$behavior <- NA
  df$behavior <- ifelse(df$nest_dist <= distance_to_nest, "nest", df$behavior)
  return(df)
}

#' Adds nest and conspecific distances
#'
#' Adds nest and conspecific distances to dataframe
#'
#' @usage AddNestConDist(df, nests_con_dist)
#'
#' @param df dataframe of locations, must have coordinate columns
#'
#' @return dataframe with "nest_dist" and "con_dist" columns
#' @export
#'
#' @details The nest_con_dist has to have structure of list[[year]][[nest_id]]
#'
AddNestConDist<- function(df,
                          nest_con_dist = nest_con_dist) {
  df <- df
  nest_con_dist <- nest_con_dist
  AddDistances <- function(df){
    df <- df
    coordinates(df) <- cbind(df$long_utm, df$lat_utm)
    proj4string(df) <- sp::CRS("+proj=utm +zone=19 +datum=NAD83")
    nest_id <- unique(df$nest_id)
    year <- unique(df$year)
    dists <- nest_con_dist[[as.character(year)]][[nest_id]]
    df$nest_dist <- round(as.vector(unlist(raster::extract(dists, df, layer=1,
      nl=1))))
    df$con_dist <- round(as.vector(unlist(raster::extract(dists, df, layer=2,
      nl=1))))
    output <- as.data.frame(df)
    return(output)
  }
  output <- ddply(df, .(year, id), AddDistances)
  output$x <- NULL
  output$y <- NULL
  return(output)
}

#' Adds nest-related columns to dataframe
#'
#' Adds nest locations, angle, and distance to each record, based on "date" and
#'   "id"
#'
#' @usage AddNestData(df, nests_study, nests_use)
#'
#' @param df dataframe of BAEA location data
#' @param nests_study .csv file of study nests. Default is:
#'   "C:/Work/R/Data/BAEA/Nests/Nests_Study.csv"
#' @param nests_use .csv file of study nests corresponding use dates. Default
#'   is: "C:/Work/R/Data/BAEA/Nests/Nests_Study_Use_Dates.csv"
#'
#' @return dataframe
#'
#' @importFrom magrittr "%>%"
#' @export
#'
AddNestData <- function(df = df,
                        nests_study = file.path("C:/Work/R/Data/BAEA/Nests",
                                        "Nests_Study.csv"),
                        nests_use = file.path("C:/Work/R/Data/BAEA/Nests",
                                      "Nests_Study_Use_Dates.csv")) {
  df <- df
  nests <- read.csv(nests_use, header=TRUE, stringsAsFactors=FALSE,
      row.names=NULL) %>%
    dplyr::left_join(.,  read.csv(nests_study, header=TRUE,
      stringsAsFactors=FALSE, row.names=NULL), by=c("name", "nest_site",
      "nest_area")) %>%
    dplyr::mutate(
      use_start_date = as.Date(as.character(use_start_date), "%Y%m%d"),
      use_end_date = as.Date(as.character(use_end_date), "%Y%m%d"),
      clutch_initiation = as.Date(as.character(clutch_initiation), "%Y%m%d"),
      breeding_end_date = as.Date(as.character(breeding_end_date), "%Y%m%d"),
      nest_long_utm = long_utm,
      nest_lat_utm = lat_utm) %>%
    dplyr::select(eagle_id, nest_site, nest_area, use_start_date, use_end_date,
      nest_long_utm, nest_lat_utm)
  nests_blank <- nests[0,]
  nests_blank[1:nrow(df),] <- NA
  df <- cbind(df, nests_blank)
    df_10000 <- df[10000,]
  for (i in 1:nrow(nests)) {
    nest <- nests[i,]
    if (is.na(nest$use_end_date)) {
      use_end_date <- Sys.Date() + 1
    } else {
      use_end_date <- nest$use_end_date
    }
    sv <- df$id == nest$eagle_id & df$date >= nest$use_start_date &
      df$date <= use_end_date
    df[sv, (ncol(df)-length(nest)+1):ncol(df)] <- nest[1,]
  }
  df <- df %>% dplyr::select(-eagle_id)
  df$nest_angle <- CalculateAngleToPoint(df$long_utm, df$lat_utm,
    df$nest_long_utm, df$nest_lat_utm)
  df <- df %>%
    dplyr::mutate(nest_dist = round(sqrt(((long_utm - nest_long_utm)^2) +
      ((lat_utm - nest_lat_utm)^2))))
  return(df)
}

#' Add perch behavior
#'
#' Finds perch behavior that meets given threshold parameters
#'
#' @usage AddPerchBehavior(df, max_speed)
#'
#' @param df dataframe
#' @param max_speed numeric, speed over which behavior is no longer considered
#'   perch.
#'
#' @return dataframe with behavior column that has "perch"
#' @export
#'
AddPerchBehavior <- function(df = df,
                             max_speed = 5) {
  df <- df
  df$behavior <- NA
  df$behavior <- ifelse(df$speed < max_speed, "perch", df$behavior)
  return(df)
}

#' Adds roost behavior to dataframe
#'
#' Finds roost arrivals and departures that meet given threshold parameters
#'
#' @usage AddRoostBehavior(df, overnight_distance_threshold,
#'   at_roost_distance_threshold, daily_location_threshold)
#'
#' @param df dataframe
#' @param default_tz used in IfElseTimedateNA/Compare functions
#' @param tz  timezone, default is "Etc/GMT+5"
#' @param overnight_distance_threshold max overnight distance
#' @param at_roost_distance_threshold  max distance away from roost
#' @param depart_timediff_max max diff between start of day and depart
#' @param arrive_timediff_max max diff between end of day and arrival
#'
#' @return dataframe with roost column that has arrive, depart, and roost
#' @export
#'
#' @details automatically makes sure that data exists for the following day,
#'   automatically checks that there are at least 7 locations in the first/last
#'   two hours of the day for departure/arrive
AddRoostBehavior <- function(df = baea,
                             default_tz = "America/New_York",
                             tz = "Etc/GMT+5",
                             overnight_distance_threshold = 100,
                             at_roost_distance_threshold = 50,
                             depart_timediff_max = 1000,
                             arrive_timediff_max = 1000){
  df$time_after_start <- as.integer(difftime(df$datetime,
    df$hr_before_sunrise, tz=tz, units = ("mins")))
  df$time_before_end <- as.integer(difftime(df$hr_after_sunset, df$datetime,
    tz=tz, units = ("mins")))
  df$two_hr_after_sunrise <- df$hr_before_sunrise + hours(2)
  df$two_hr_before_sunset <- df$hr_after_sunset - hours(2)
  df$sunrise_window_loc <- df$datetime <= df$two_hr_after_sunrise
  df$sunset_window_loc <- df$datetime >= df$two_hr_before_sunset
  sumstats <- ddply(df, .(id, date), summarize,
    date = as.Date(unique(date)), total_loc = length(deploy_seq),
    am_loc = sum(sunrise_window_loc, na.rm=TRUE),
    pm_loc = sum(sunset_window_loc, na.rm=TRUE))
  nextday <- function(data = data){
    out <- sapply(2:nrow(data),function(i){data$date[i] - data$date[i-1]})
    next_day <- c(out, NA)
    next_day
  }
  nextamloc <- function(data = data){
    out <- sapply(1:nrow(data),function(i) { data$am_loc[i+1] })
    next_am_loc <- c(out)
    next_am_loc
  }
  list <- by(sumstats, sumstats$id, function(x) nextday(x))  # makes list
  sumstats$nextdayGPS <- unlist(list)
  list <- by(sumstats, sumstats$id, function(x) nextamloc(x))  # makes list
  sumstats$next_am_loc <- unlist(list)
  # At this point, sumstats has: "id", "date", "total_locs", "nextdayGPS",
  # "am_loc", "pm_loc", and "next_am_loc"
  # This has all the data needed to cull by nextDayGPS, and am/pm locations
  df<-merge(df, sumstats, by = c("id", "date"), all.x = TRUE)
  last_threshold<- function(data, threshold=overnight_distance_threshold){
    subset (data, last == "Last" & step_length <= threshold & nextdayGPS == 1 &
    pm_loc >= 7 & next_am_loc >= 7, select=c("id", "date"))
  }
  last_roost_confirmed <- last_threshold(data = df,
    threshold=overnight_distance_threshold)
  # This culls by overnight segment length, if next day has GPS locations, if
  # there are at least 7 locations within an hour of either side of sunset
  # and if there are at least 7 locations within an hour of either side of
  # sunrise on the following morning.
  row.names(last_roost_confirmed) <- NULL  # housekeeping
  first_roost_confirmed <- adply(last_roost_confirmed, 1, transform,
    date = date+1)  # the mornings after "last_roost_confirmed" dates
  roost_arrival_filtered <- join (df, last_roost_confirmed, type="inner")
    # to "confirm" arrival based on overnight distance
  roost_departure_filtered <- join (df, first_roost_confirmed,
    type = "inner")  # to "confirm" departure based on overnight distance
  threshold_dist <- function(x, threshold=at_roost_distance_threshold) {
    x>threshold
  }
  # This function sets the distance a location can still be considered at roost
  # based on its distance to the last and first locations
  depart <- ddply(roost_departure_filtered, .(date, id), function(x)
    x[(Position(threshold_dist, x$dist_first, right = FALSE, nomatch = NULL)
       - 1), c("id","date", "datetime")])
  arrive <- ddply(roost_arrival_filtered, .(date, id), function(x)
    x[(Position(threshold_dist, x$dist_last, right = TRUE, nomatch = NULL) + 1),
      c("id", "date", "datetime")])
  arrive$arr_dist_threshold_datetime <- arrive$datetime
  depart$dep_dist_threshold_datetime <- depart$datetime
  arrive$datetime <- NULL
  depart$datetime <- NULL
  threshold_time_depart <- function(x, threshold=depart_timediff_max) {
    x>threshold
  }
  threshold_time_arrive <- function(x, threshold=arrive_timediff_max) {
    x>threshold
  }
  depart_max <- ddply(roost_departure_filtered, .(date, id),
    function(x) x[nrow(x), c("id","date", "datetime")])
  arrive_max <- ddply(roost_arrival_filtered, .(date, id),
    function(x) x[1, c("id","date", "datetime")])
  depart_max$max_datetime <- depart_max$datetime
  arrive_max$max_datetime <- arrive_max$datetime
  depart_max$datetime <- NULL
  arrive_max$datetime <- NULL
  # max_datetime is the first/last record in the day
  depart_threshold <- ddply(roost_departure_filtered, .(date, id), function(x)
    x[(Position(threshold_time_depart, x$time_after_start, right = FALSE,
    nomatch = NULL) - 1), c("id","date", "datetime")])
  arrive_threshold <- ddply(roost_arrival_filtered, .(date, id), function(x)
    x[(Position(threshold_time_arrive, x$time_before_end, right = TRUE,
    nomatch = NULL) + 1), c("id","date", "datetime")])
  depart_threshold$threshold_datetime <- depart_threshold$datetime
  arrive_threshold$threshold_datetime <- arrive_threshold$datetime
  depart_threshold$datetime <- NULL
  arrive_threshold$datetime <- NULL
  # threshold_dateime is the first record within the threshold
  depart_max <- merge(depart_max, depart_threshold, by = c("date", "id"),
    all.x= TRUE)
  arrive_max <- merge(arrive_max, arrive_threshold, by = c("date", "id"),
    all.x= TRUE)
  rm(arrive_threshold, depart_threshold)
  depart_max <- IfElseTimedateNA(df=depart_max, col1="threshold_datetime",
    col2="max_datetime", result="max_threshold", default_tz=default_tz, tz=tz)
  arrive_max <- IfElseTimedateNA(df=arrive_max, col1="threshold_datetime",
    col2="max_datetime", result="max_threshold", default_tz=default_tz, tz=tz)
  depart <- merge(depart, depart_max, by = c("date", "id"), all.x= TRUE)
  arrive <- merge(arrive, arrive_max, by = c("date", "id"), all.x= TRUE)
  depart <- IfElseTimedateCompare(df=depart, col1="dep_dist_threshold_datetime",
    sign="<", col2="max_threshold", result="dep_datetime",
    default_tz=default_tz, tz=tz)
  arrive <- IfElseTimedateCompare(df=arrive, col1="max_threshold", sign=">",
    col2="arr_dist_threshold_datetime", result="arr_datetime",
    default_tz=default_tz, tz=tz)
  roost_times <- merge(depart, arrive, by = c("date", "id"), all = TRUE)
  roost_times <- subset(roost_times, select=c("date", "id", "dep_datetime",
    "dep_dist_threshold_datetime", "arr_datetime",
    "arr_dist_threshold_datetime"))
  df <- merge(df, roost_times, by = c("date", "id"), all= TRUE)
  df <- df[with(df, order(id, date, datetime)), ]
  df$loaf <- NA
  df$loaf <- ifelse(df$datetime <= df$dep_dist_threshold_datetime, "loaf", NA)
  df$loaf <- ifelse(is.na(df$loaf) & df$datetime >=
    df$arr_dist_threshold_datetime, "loaf", df$loaf)
  df$depart <- ifelse(df$datetime == df$dep_datetime, "depart", NA)
  df$arrive <- ifelse(df$datetime == df$arr_datetime, "arrive", NA)
  df$depart <- ifelse(is.na(df$depart) & !is.na(df$dep_datetime) &
                         df$datetime < df$dep_datetime, "roost", df$depart)
  df$arrive <- ifelse(is.na(df$arrive) & !is.na(df$arr_datetime) &
                         df$datetime > df$arr_datetime, "roost", df$arrive)
  df$roost <- ifelse(is.na(df$depart),df$arrive,df$depart)
  df$roost_loaf <- ifelse(!is.na(df$roost), df$roost,
    df$loaf)
  if (!("behavior" %in% colnames(df))) df$behavior<-NA
  df$behavior <- ifelse(!is.na(df$roost_loaf) & is.na(df$behavior),
    df$roost_loaf, df$behavior)
  drops <- c("arrive", "depart", "total_loc", "nextdayGPS", "dep_datetime",
    "two_hr_after_sunrise", "two_hr_before_sunset", "sunrise_window_loc",
    "arr_datetime", "sunset_window_loc", "am_loc", "pm_loc", "next_am_loc",
    "roost", "loaf", "roost_loaf" ,"dep_dist_threshold_datetime",
    "arr_dist_threshold_datetime")
  df<-df[ ,!(names(df) %in% drops)]
  row.names(df) <- NULL
  return(df)
}

#' Adds time step proportions
#'
#' Adds time_steps, day_min, time_after_start, and time_proportion data
#'
#' @param df dataframe
#' @param by grouping column. Default is "id"
#' @param time_step time step length, based on lubriate times. Default is
#'   "15 min".
#' @param tz timezone for dataset. Default is "Etc/GMT+5"
#'
#' @return dataframe with time_steps, day_min, time_after_start, and
#'   time_proportion data
#' @export
#'
#' @details  need to run AddSolarTimes() prior to running this function
#'
AddTimeStepProportion <- function(df = df,
                                  by = "id",
                                  time_step = "15 min",
                                  tz = "Etc/GMT+5") {
  df <- df
  df$by <- df[,by]
  DailyTimeStepCount <- function (df=df){
    days <- subset(df, select=c("id", "date", "hr_before_sunrise",
      "hr_after_sunset"))
    days <- plyr::ddply(days, plyr::.(date), function(x) x[1, ]) # first records
    days$time_steps <- NA
    days$day_min <- NA
    for (i in 1:nrow(days)) {
      day <- days[i,]
      days[i,"time_steps"] <- length(seq(day[,"hr_before_sunrise"],
        day[,"hr_after_sunset"], time_step))
      days[i,"day_min"] <- as.integer(difftime(day[,"hr_after_sunset"],
        day[,"hr_before_sunrise"],tz=tz, units = ("mins")))
    }
  return(days)
  }
  df2 <- plyr::ddply(df, plyr::.(by), DailyTimeStepCount)
  df2 <- merge(df, df2, all.x=TRUE)
  df2$time_after_start <- as.integer(difftime(df2$datetime,
    df2$hr_before_sunrise, tz=tz, units = ("mins")))
  df2$time_proportion <- df2$time_after_start/df2$day_min
  df2$by <- NULL
  return(df2)
}

#' Analyses overlap of homeranges in paired nest locations
#'
#' Analyzes pair locations to determine mid-way point betweeen nests and
#'   proportion of overlapping locations
#'
#' @usage AnalyzePairLocation(pair, nests, points, mid_length, zoom)
#'
#' @param pair pair to be analyzed, first nest is analyzed
#' @param nests nest locations
#' @param points baea locations
#' @param mid_length length of mid-line.
#' @param zoom  integer (2-22) of google map zoom, default = 10.
#'
#' @return Multiple output: writes .csv, plots, and dataframe
#'
#' @importFrom magrittr "%>%"
#' @export
#'
AnalyzePairLocations <- function(pair,
                                 nests,
                                 points,
                                 mid_length,
                                 zoom = 12){
  crs_utm <- "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  crs_proj <- "+proj=longlat +datum=WGS84"
  nests_pair <- nests[nests$name %in% pair,]
  nests_pair_df <- as.data.frame(nests_pair)
  pair_dist <- round(sqrt(sum((c(nests_pair$long_utm[1], nests_pair$lat_utm[1])-
      c(nests_pair$long_utm[2], nests_pair$lat_utm[2]))^2)))/1000
  points1 <- points[points$id == pair[1], ]
  mid_pt_df <- cbind.data.frame(long_utm = (sum(nests_pair$long_utm))/2,
    lat_utm = (sum(nests_pair$lat_utm))/2)
  mid_pt <- SpatialPointsDataFrame(mid_pt_df[c("long_utm", "lat_utm")],
    bbox=NULL, data=mid_pt_df, proj4string=CRS(crs_utm))
  angle1 <- ConvertAngle(CalculateAngleToPoint(nests_pair$long_utm[1],
    nests_pair$lat_utm[1], nests_pair$long_utm[2], nests_pair$lat_utm[2]) +
    (pi/2))
  angle2 <- ConvertAngle(CalculateAngleToPoint(nests_pair$long_utm[1],
    nests_pair$lat_utm[1], nests_pair$long_utm[2], nests_pair$lat_utm[2]) -
    (pi/2))
  mid_ends_df <- data.frame(long_utm=NA, lat_utm=NA)
  mid_ends_df[1,]<- CoordinatesFromAngleLength(mid_pt_df[1,], angle1,mid_length)
  mid_ends_df[2,]<- CoordinatesFromAngleLength(mid_pt_df[1,], angle2,mid_length)
  mid_ends <- SpatialPointsDataFrame(mid_ends_df[c("long_utm", "lat_utm")],
    bbox=NULL, data=mid_ends_df, proj4string=CRS(crs_utm))
  mid_line <- CreateSpatialLines(df=mid_ends_df, long="long_utm", lat="lat_utm",
    crs=crs_utm)
  nests_pair_df <- as.data.frame(nests_pair) %>% dplyr::select(long_utm,lat_utm)
  poly_df <- rbind(nests_pair_df[1,], mid_ends_df[1,], nests_pair_df[2,],
    mid_ends_df[2,])
  poly_sp <- CreateSpatialPolygons(poly_df, long="long_utm", lat="lat_utm",
    crs=crs_utm)
  points1_inside <- as.data.frame(points1[poly_sp, "home_dist"])
  points1_inside_count <- nrow(points1_inside)
  points1_inside$home_nest <- pair[1]
  points1_inside$con_nest <- pair[2]
  points1_inside$pair_dist <- pair_dist
  points1_inside$mid_length <- mid_length
  points1_inside <- points1_inside %>% dplyr::select(home_nest, pair_dist,
    everything(), con_nest)
  out_folder <- file.path(getwd(), "Points_Inside", mid_length)
  ifelse(!dir.exists(out_folder), dir.create(out_folder), FALSE)
  write.csv(points1_inside, file=file.path(out_folder, paste0(pair[1], "_",
    pair[2], ".csv")), row.names = FALSE)
  poly_sp_df <- tidy(poly_sp)
  centroids <- coordinates(poly_sp)
  x <- centroids[,1]
  y <- centroids[,2]
  poly_spdf <- SpatialPolygonsDataFrame(poly_sp, data=data.frame(x=x, y=y,
    row.names=row.names(poly_sp)))
  rgdal::writeOGR(obj=poly_spdf,  dsn = file.path(getwd(), "Polys"),
    layer=paste0(pair[1], "_", pair[2]), driver = "ESRI Shapefile",
    overwrite_layer = TRUE)

  out_folder <- file.path(image_output, "Pair Maps", mid_length)
  ifelse(!dir.exists(out_folder), dir.create(out_folder), FALSE)

  out_folder <- file.path(image_output, "Pair Locations", mid_length)
  ifelse(!dir.exists(out_folder), dir.create(out_folder), FALSE)

  # Location Graphs
  gg <- ggplot() + theme_no_legend +
    geom_polygon(data=poly_sp_df, aes(x=long, y=lat), fill="dodgerblue1",
      colour="black", size=1) +
    geom_path(data=as.data.frame(coordinates(mid_line)), aes(x=long_utm,
      y=lat_utm), colour="grey80", size=1, linetype=2)
  x_diff <- diff(range(poly_sp_df$long))
  y_diff <- diff(range(poly_sp_df$lat))
  axis_max <- which.max(c(x_diff, y_diff))
  if (axis_max == 1) {
    y_mean <- mean(range(poly_sp_df$lat))
    y_add <- diff(range(poly_sp_df$long)/2)
    gg <- gg + coord_fixed(ratio=1, ylim=c(y_mean - y_add, y_mean + y_add),
      xlim = range(poly_sp_df$long))
  } else {
    x_mean <- mean(range(poly_sp_df$long))
    x_add <- diff(range(poly_sp_df$lat)/2)
    gg <- gg + coord_fixed(ratio=1, xlim=c(x_mean - x_add, x_mean + x_add),
      ylim = range(poly_sp_df$lat))
  }
  gg <- gg + geom_point(data=as.data.frame(points1), aes(x=long_utm,y=lat_utm),
      colour="red", size=1) +
    geom_point(data=points1_inside, aes(x=long_utm,y=lat_utm),
      colour="yellow", size=1.25) +
    ggtitle(paste0(pair[1], " - ", pair[2], " (n=", points1_inside_count,
      ")"))
  plot(gg)
  SaveGGPlot(paste0(pair[1], "_", pair[2], ".jpeg"), file.path(image_output,
    "Pair Locations", mid_length))

  # GG Maps
  points1_inside_sp <- points1[poly_sp, "home_dist"]
  points1_inside_wgs84 <- as.data.frame(spTransform(points1_inside_sp,
    CRS(crs_proj)))
  points1_wgs84 <- as.data.frame(spTransform(points1, CRS(crs_proj)))
  mid_pt_wgs84 <- as.data.frame(spTransform(mid_pt, CRS(crs_proj)))
  mid_ends_wgs84 <- spTransform(mid_ends, CRS(crs_proj))
  poly_wgs84 <- tidy(spTransform(poly_sp, CRS(crs_proj)))
  points1_inside_count <- nrow(points1_inside_wgs84)
  ptspermm <- 2.83464567
  title_text <- paste0(pair[1], " - ", pair[2], " (n=", points1_inside_count,
    ")")
  nest_map <- get_map(location = c(lon=mid_pt_wgs84[1,3],
    lat=mid_pt_wgs84[1,4]), color="color", source="google",
    maptype="hybrid", zoom=zoom)
  nest_ggmap <- ggmap(nest_map, extent="device", ylab="Latitude",
    xlab="Longitude", legend="right")
#  Sys.sleep(2)
  bb <- attr(nest_map, "bb")
#  bb_diag <- bb$ll.lat - bb$ur.lat
  scalebar.length = 5
  sbar <- data.frame(long_start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                     long_end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                     lat_start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                     lat_end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                     bb_diag = bb$ur.lat - bb$ll.lat,
                     scalebar.length = scalebar.length)
  map_title <- cbind.data.frame(name=title_text, x = (bb$ll.lon + bb$ur.lon)/2,
    y = bb$ur.lat - 0.01)
  sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$long_start,
    sbar$lat_start), c(sbar$long_end,sbar$lat_end))
  sbar$long_end <- sbar$long_start +
    ((sbar$long_end-sbar$long_start)/sbar$distance) * scalebar.length*1000
  gg <-
    nest_ggmap +
    geom_polygon(data=poly_wgs84, aes(x=long, y=lat), alpha=0.8,
      fill="dodgerblue1", color="black", size=0.2) +
    geom_path(data=as.data.frame(coordinates(mid_ends_wgs84)), aes(x=long_utm,
      y=lat_utm), colour="grey80", size=1, linetype=2) +
    geom_point(data=points1_wgs84, aes(x=long, y=lat), shape=19, alpha=.9,
      color="red", size=.75, stroke=1, na.rm=TRUE) +
    geom_point(data=points1_inside_wgs84, aes(x=long_utm,y=lat_utm), shape=19,
      alpha=.75, color="yellow", size=2, stroke=1, na.rm=TRUE) +
   geom_segment(data = sbar, aes(x=long_start, xend=long_end, y=lat_start,
      yend=lat_end), size=1.25, arrow=arrow(angle=90, length=unit(0.2, "cm"),
      ends="both", type="open"), color="white") +
    geom_text(data = sbar, aes(x = (long_start + long_end)/2, y = lat_start +
      0.025*(bb_diag), label=paste(format(scalebar.length),'km')),
      hjust=0.5, vjust=.7, size=18/ptspermm, color="white")  +
    geom_text(data = map_title, aes(x=x, y=y, label=name), size=24/ptspermm,
      color="white") +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  print(gg)
  SaveGGPlot(paste0(pair[1], "_", pair[2], ".jpeg"), file.path(image_output,
    "Pair Maps", mid_length), bg = "black")
  return(points1_inside)
}

#' Converts nest alphanumeric id to nest numeric id
#'
#' Converts alphanumeric "nest_id" column to a numeric "nest_id_num" column,
#'   useful for creating RasterLayers associated with nests
#'
#' @usage ConvertNestIdToNum(df)
#'
#' @param df input dataframe with "nest_id" column
#'
#' @return A dataframe with a numberic "nest_id_num" column
#' @export
#'
#' @details  If a "nest_id_num" column doesn't exist, the function
#'   automatically creates one.
#'
ConvertNestIdToNum <- function(df){
  df <- df
  if(!"nest_id" %in% colnames(df)) {
    df$nest_id <- df$nest_site
  }
  if(!"nest_id_num" %in% colnames(df)) {
    df$nest_id_num <- NA
  }
  nest_id <- df$nest_id
  territory_number <- sapply(strsplit(nest_id, "[A-Z]"), "[", 1)
  nest_letter <- stringr::str_extract(nest_id, "[A-Z]")
  nest_number <- sapply(strsplit(nest_id, "[A-Z]"), "[", 2)
  nest_number[is.na(nest_number)] <- ""
  letters <- LETTERS[1:26]  # had only gone to "k" in 2012
  numbers <- formatC(1:26, width = 2, format = "d", flag = "0")
  for (i in 1:length(letters)){
    nest_letter <- gsub(letters[i], numbers[i], nest_letter)
  }
  df$nest_id_num <- as.numeric(sprintf("%s%s%s", territory_number, nest_letter,
    nest_number))
  return(df)
}

#' Converts nest numeric id to nest alphanumeric id
#'
#' Converts a numeric "nest_id_num" column to an alphanumeric "nest_id" column
#'
#' @usage ConvertNestIdToNum(df)
#'
#' @param df input dataframe with "nest_id_num" column
#'
#' @return A dataframe with a numberic "nest_id" column
#' @export
#'
#' @details If a "nest_id" column doesn't exist, the function automatically
#'   creates one
ConvertNestNumToId <- function(df){
  df <- df
  if(!"nest_id" %in% colnames(df)) {
    df$nest_id <- NA
  }
  nest_id_num <- df$nest_id_num
  letters <- LETTERS[1:26]  # had only gone to "k" in 2012
  numbers <- formatC(1:26, width = 2, format = "d", flag = "0")
  for (i in 1:length(nest_id_num)){
    if (nchar(nest_id_num[i]) == 3){
      territory_number<- sprintf("%03d", as.numeric(substr(nest_id_num[i],1,1)))
      nest_letter <- substr(nest_id_num[i], 2, 3)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, sep="")
    }
    if (nchar(nest_id_num[i]) == 4){
      territory_number<- sprintf("%03d", as.numeric(substr(nest_id_num[i],1,2)))
      nest_letter <- substr(nest_id_num[i], 3, 4)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, sep="")
    }
    if (nchar(nest_id_num[i]) == 5){
      territory_number <- substr(nest_id_num[i], 1, 3)
      nest_letter <- substr(nest_id_num[i], 4, 5)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, sep="")
    }
    if (nchar(nest_id_num[i]) == 7){
      territory_number <- substr(nest_id_num[i], 1, 3)
      nest_letter <- substr(nest_id_num[i], 4, 5)
      nest_number <- substr(nest_id_num[i], 6, 7)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, nest_number,
        sep="")
    }
    if (nchar(nest_id_num[i]) == 8){
      territory_number <- substr(nest_id_num[i], 1, 4)
      nest_letter <- substr(nest_id_num[i], 5, 6)
      nest_number <- substr(nest_id_num[i], 7, 8)
      for (j in 1:length(letters)){
        nest_letter <- gsub(numbers[j], letters[j], nest_letter)
      }
      df$nest_id[i] <- paste(territory_number, nest_letter, nest_number,
        sep="")
    }
  }
  return(df)
}

#' Convert step data coordinates to WGS84
#'
#' Converts "x", "y" columns to coordinates to lat/long WGS84 (for Google)
#'
#' @usage ConvertStepDataCoordinates(df)
#'
#' @param df input dataframe with "x", "y" columns
#' @param crs string of projection for "x", "y" columns, default is for Maine
#'   BAEA GPS data (UTM Zone 19N).
#'
#' @return  A dataframe with added "long", "lat" columns
#' @export
#'
#' @details
#'
ConvertStepDataCoordinates <- function(df,
                                       crs = "+proj=utm +zone=19 ellps=WGS84"){

  df <- df
  sp::coordinates(df) <- c("x", "y")
  sp::proj4string(df) <- sp::CRS(crs)
  res <- sp::spTransform(df, CRS("+proj=longlat +datum=WGS84"))
  long_lat <- sp::coordinates(res)
  colnames(long_lat) <- c("long", "lat")
  output <- cbind.data.frame(df, long_lat)
  output$optional <- NULL
  return(output)
}

#' Creates colors based on any factor in the dataframe
#'
#' Creates and/or displays dataframe of "by" variable and associated colors#'
#'
#' @usage CreateColorsByAny(by, df, r_pal, b_pal, output, plot, ...)
#'
#' @param by variable to determine colors. Specific outcomes for: "behavior",
#'   "id", or "sex"
#' @param df dataframe with "by" variable - only required if "by" is not
#'   "behavior", "id", or "sex"
#' @param pal dataframe with "by" variable - only required if "by" is not
#'   "behavior", "id", or "sex"
#' @param r_pal color palette for CreateVarsColors(), default is NULL
#' @param b_pal color palette from RColorBrewer for CreateVarsColors(), default
#'   is "Set1".
#' @param output color palette from RColorBrewer for CreateVarsColors(),
#'   default is "Set1"
#' @param plot logical, whether or not to display names and colors, default is
#'   FALSE
#' @param ... additional parameters for CreateColorsByMetadata()
#'
#' @return df with "by" variable and hexidecimal colors
#' @export
#'
#' @details Used in several other functions. Color palettes determined
#'   automatically for "behavior", "id", and "sex". Others set by pal or b_pal.
#'   Requires 'kml' associated functions.
CreateColorsByAny <- function (by,
                               df,
                               pal = NULL,
                               r_pal = NULL,
                               b_pal = "Set1",
                               output = TRUE,
                               plot = FALSE,
                               ...) {
  if (!is.null(by)){
  if (by == "behavior" || by == "id" || by == "sex") {
    if (by == "behavior") by_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/Models/Behavior_Colors.csv", metadata_id="behavior")
    if (by == "id") by_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/GPS_Deployments.csv", metadata_id="deploy_location")
    if (by == "sex") by_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/Models/Behavior_Colors.csv", metadata_id="sex")
  } else {
    by_colors <- CreateColorsByVar(df=df, by=by, pal=pal, r_pal = r_pal, b_pal =
      b_pal)
  }
  } else {
    by_colors <- CreateColorsByVar(df=df, by=by, pal=pal, r_pal = r_pal, b_pal =
      b_pal)
  }
  if (plot == TRUE) PlotColorPie(by_colors)
  if (output == TRUE) return(by_colors)
}

#' Create homerange kernels for the BAEA nests
#'
#' Creates RasterLayers of homerange kernels based on home centroids
#'
#' @usage CreateHomeRangeKernelsBAEA(df_all, df_home, base, max_r,
#'   home_inflection, home_scale, avoid_inflection, avoid_scale, output_dir,
#'   id, write_distance, write_homerange)
#'
#' @param df_all dataframe of all homerange centroids
#' @param df_home dataframe of homerange centroids to calculate homerange
#'   kernels (these must be a subset of the df_all dataframe), Default is to
#'   use 'df_all' dataframe
#' @param base base Raster that sets the projection, extent, and dimensions of
#'   the study area
#' @param max_r maximum radius to calculate the homerange raster from each
#'   df_home centroid
#' @param home_inflection inflection point of the Logistic function that
#'   governs the home kernel
#' @param home_scale scale parameter of the Logistic function that governs the
#'   scale parameter
#' @param avoid_inflection inflection point of the Logistic function that
#'   governs the conspecific avoidance kernel
#' @param avoid_scale scale parameter of the Logistic function that governs the
#'   conspecific avoidance kernel
#' @param min_prob numeric, minimum probablity threshold; all probabilities
#'   below this value will be rounded to zero. Default is 0.
#' @param id column name of df_home that identifies the homerange. Default is
#'   NULL, which sets the names to df_home row number.
#' @param name column name of df_home that identifies the name of the homerange
#'   and is used to write the file name of the .tif Default is 'id', which sets
#'  the names to df_home row number if 'id' is NULL.
#' @param output_dir directory for output files (distance, homerange)
#' @param write_distance logical, write distance raster to file. Default is
#'   FALSE.
#' @param write_homerange logical, write home range raster to file. Default is
#'   FALSE.
#'
#' @return A list containing homerange kernel Rasters for all the df_home
#'   centroids.
#' @export
#'
#'
CreateHomeRangeKernelsBAEA <- function(df_all,
                                       df_home = df_all,
                                       base,
                                       max_r,
                                       home_inflection,
                                       home_scale,
                                       avoid_inflection,
                                       avoid_scale,
                                       min_prob,
                                       id = NULL,
                                       name = id,
                                       output_dir,
                                       write_distance = FALSE,
                                       write_homerange = FALSE) {
  cellsize <- raster::res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- raster::xmin(base)
  ymin <- raster::ymin(base)
  df_all_sp <- sp::SpatialPointsDataFrame(df_all[,c("x","y")], df_all,
    proj4string=crs(base))
  homerange_ids <- df_home[,id]
  total <- length(homerange_ids)
  homerange_kernels <- as.list(setNames(rep(NA,length(homerange_ids)),
    homerange_ids), homerange_ids)
  for (j in 1:nrow(df_home)) {
    home <- df_home[j,]
    ifelse(is.null(name), home_name <- j, home_name <- home[,name])
    writeLines(noquote(paste0("Calculating homerange for: ", home_name,
      " (", j, " of ", total, ").")))
    home_sp <- sp::SpatialPointsDataFrame(home[,c("x","y")], home,
      proj4string=crs(base))
    xy <- CenterXYInCell(home_sp@coords[,"x"], home_sp@coords[,"y"],
      xmin, ymin, cellsize)
    cell_extent <- raster::extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
      (cellsize/2), xy[2]+(cellsize/2))
    cell <- raster::setValues(raster::raster(cell_extent,
      crs=raster::projection(base), res=cellsize), j)
    home_ext <- raster::extend(cell, c(max_r_cells, max_r_cells), value=NA)
    home_dist <- raster::distanceFromPoints(home_ext, home[,c("x","y")])
    home_kern <- raster::calc(home_dist, fun = function(x){(1/(exp((-(x -
      home_inflection)) / home_scale) + 1))})
    base_crop <- raster::crop(base, extent(home_ext))
    global_dist_crop <- raster::distanceFromPoints(base_crop, df_all_sp)
    if (write_distance == TRUE) {
      filename <- paste0(output_dir,"/ConDist_", home_name, ".tif")
      global_dist_crop <- raster::distanceFromPoints(base_crop, df_all_sp)
      raster::writeRaster(global_dist_crop, filename=filename, format="GTiff",
        overwrite=TRUE)
      writeLines(noquote(paste("Writing:", filename)))
    } else {
      global_dist_crop <- raster::distanceFromPoints(base_crop, df_all_sp)
    }
    cent_dist <- raster::overlay(home_dist, global_dist_crop, fun=function (x,y)
      {ifelse(x != y, NA, x)})
    cent_bounds <- raster::boundaries(cent_dist)
    cent_bounds <- raster::subs(cent_bounds, data.frame(from=c(0,1),
      to=c(NA,1)))
    edge_dist_abs <- raster::distance(cent_bounds)
    edge_dist <- raster::overlay(cent_dist, edge_dist_abs, fun=function(x,y)
      {ifelse(!is.na(x), y*-1, y*1)})
    avoid_kern <- raster::calc(edge_dist, fun = function(x){(1/(exp((-(x -
      avoid_inflection)) / avoid_scale) + 1))})
    homerange_kern <- raster::overlay(avoid_kern, home_kern, fun=function(x,y){
      p <- y*x
      })
    homerange_kern[homerange_kern < min_prob] <- 0
    if (write_homerange == TRUE) {
      filename <- paste0(output_dir,"/HomeRange_",home_name, ".tif")
      writeLines(noquote(paste0("Writing: ", filename)))
      raster::writeRaster(homerange_kern, filename=filename, format="GTiff",
        overwrite=TRUE)
    }
    homerange_kernels[[j]] <- homerange_kern
    names(homerange_kernels[[j]]) <- df_home[j, id]
  }
  return(homerange_kernels)
}

#' Create homerange kernels based on interpolating empiricial x,y location data
#'
#' Creates RasterLayers of homerange kernels based on home centroids
#'
#' @usage CreateHomeRangeKernelsInterpolate(df_all, df_home, base, max_r, id,
#'   name, min_prob, output_dir, write_homerange)
#'
#' @param df_all dataframe of all homerange centroids
#' @param df_home dataframe of homerange centroids to calculate homerange
#'   kernels (these must be a subset of the df_all dataframe), Default is to
#'   use 'df_all' dataframe
#' @param base base Raster that sets the projection, extent, and dimensions of
#'   the study area
#' @param max_r maximum radius to calculate the homerange raster from each
#'   df_home centroid
#' @param id column name of df_home that identifies the homerange.
#' @param name column name of df_home that identifies the name of the homerange
#'   and is used to write the name of the .tif file.  Default is 'id', which
#'   sets the names to df_home row number if 'id' is NULL.
#' @param min_prob numeric, minimum probablity threshold; all probabilities
#'   below this value will be rounded to zero.  Default is 0.
#' @param output_dir directory for output files (distance, homerange)
#' @param write_distance logical, write distance range raster to file. Default
#'   is FALSE
#' @param write_homerange logical, write home range raster to file. Default is
#'   FALSE
#'
#' @return A list containing homerange kernel Rasters for all the df_home
#'   centroids
#' @export
#'
#' @details Used inside AddHomeRangeKernelsBAEA()
#'
CreateHomeRangeKernelsInterpolate <- function(df_all,
                                              df_home = df_all,
                                              base,
                                              max_r,
                                              id,
                                              name = NULL,
                                              min_prob = 0,
                                              output_dir,
                                              write_distance = FALSE,
                                              write_homerange = FALSE) {
  cellsize <- raster::res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- raster::xmin(base)
  ymin <- raster::ymin(base)
  df_all_sp <- SpatialPointsDataFrame(df_all[,c("x","y")], df_all,
    proj4string=raster::crs(base))
  homerange_ids <- df_home[,id]
  total <- length(homerange_ids)
  homerange_kernels <- as.list(setNames(rep(NA,length(homerange_ids)),
    homerange_ids), homerange_ids)
  for (j in 1:nrow(df_home)) {
    home <- df_home[j,]
    ifelse(is.null(name), home_name <- j, home_name <- home[,name])
    writeLines(noquote(paste0("Calculating homerange for: ", home_name,
      " (", j, " of ", total, ").")))
    home_sp <- sp::SpatialPointsDataFrame(home[,c("x","y")], home,
      proj4string=raster::crs(base))
    xy <- CenterXYInCell(home_sp@coords[,"x"], home_sp@coords[,"y"],
      xmin, ymin, cellsize)
    cell_extent <- raster::extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
      (cellsize/2), xy[2]+(cellsize/2))
    cell <- setValues(raster::raster(cell_extent, crs=raster::projection(base),
      res=cellsize),j)
    home_ext <- raster::extend(cell, c(max_r_cells, max_r_cells), value=NA)
    ### NEW SECTION STARTS HERE ###
    df_all$r_value <- 0
    df_all[df_all[,id] == home[1,id], "r_value"] <- 1
    all <- CreateRasterFromPointsAndBase(df_all, "r_value", "x", "y",
      base=home_ext)
    xy <- data.frame(xyFromCell(all, 1:ncell(all)))
    xy_values <- raster::getValues(all)
    # thin plate spline model
    tps_model <- fields::Tps(xy, xy_values)
    # use model to predict values at all locations
    homerange_kern <- raster::interpolate(home_ext, tps_model)
    ### NEW SECTION ENDS HERE ###
    homerange_kern[homerange_kern < min_prob] <- 0
    if (write_homerange == TRUE) {
      filename <- paste0(output_dir,"/HomeRange_",home_name, ".tif")
      writeLines(noquote(paste0("Writing: ", filename)))
      writeRaster(homerange_kern, filename=filename, format="GTiff",
        overwrite=TRUE)
      ExportKMLRasterOverlay(homerange_kern, alpha =.8, outfile=paste0(name,
        "_Interpolate"), output_dir=getwd())
    }
    homerange_kernels[[j]] <- homerange_kern
    names(homerange_kernels[[j]]) <- df_home[j, name]
  }
  return(homerange_kernels)
}

#' Create home range kernels using a Pareto distribution
#'
#' Creates RasterLayers of homerange kernels based on home centroids
#'
#' @usage CreateHomeRangeKernelsParetoGamma(df_all, df_home, base, max_r,
#'   home_shape, home_scale, min_prob, id, name, output_dir, write_homerange)
#'
#' @param df_all dataframe of all homerange centroids
#' @param df_home dataframe of homerange centroids to calculate homerange
#'   kernels (these must be a subset of the df_all dataframe), Default is to
#'   use 'df_all' dataframe.
#' @param base base Raster that sets the projection, extent, and dimensions of
#'   the study area
#' @param max_r maximum radius to calculate the homerange raster from each
#'   df_home centroid
#' @param home_shape inflection point of the Logistic function that governs the
#'   home kernel
#' @param home_scale scale parameter of the Logistic function that governs the
#'   scale parameter
#' @param id column name of df_home that identifies the homerange. Default is
#'   NULL, which sets the names to df_home row number.
#' @param name column name of df_home that identifies the name of the homerange
#'   and is used to write the name of the .tif file. Default is 'id', which
#'   sets the names to df_home row number if'id' is NULL.
#' @param min_prob numeric, minimum probablity threshold; all probabilities
#'   below this value will be rounded to zero.  Default is 0.
#' @param output_dir directory for output files (distance, homerange)
#' @param write_homerange logical, write home range raster to file. Default is
#'   FALSE
#'
#' @return A list containing homerange kernel Raster for all the df_home
#'   centroids
#' @export
#'
CreateHomeRangeKernelsPareto <- function(df_all,
                                         df_home = df_all,
                                         base = base,
                                         max_r,
                                         home_shape,
                                         home_scale,
                                         id = "nest_site",
                                         name ="nest_id",
                                         min_prob = 0,
                                         output_dir = getwd(),
                                         write_homerange = TRUE){
  source('C:/Work/R/Functions/sim/move.R')
  cellsize <- res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  df_sp <- SpatialPointsDataFrame(df_all[,c("x","y")], df_all,
    proj4string=crs(base))
#  base <- crop(base, df_sp, snap='out')
#  base <- extend(base, c(max_r_cells/2, max_r_cells/2), value=11)
  homerange_ids <- df_home$nest_id
  total <- length(homerange_ids)
  homerange_kernels <- as.list(setNames(rep(NA,length(homerange_ids)),
    homerange_ids), homerange_ids)
  for (j in 1:nrow(df_home)) {
    writeLines(noquote(paste("Calculating homerange", j, "of", total)))
    home <- df_home[j,]
    home_sp <- SpatialPointsDataFrame(home[,c("x","y")], home,
      proj4string=crs(base))
    xy <- CenterXYInCell(home_sp@coords[,"x"], home_sp@coords[,"y"],
      xmin, ymin, cellsize)
    cell_extent <- extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
      (cellsize/2), xy[2]+(cellsize/2))
    cell <- setValues(raster(cell_extent, crs=projection(base), res=cellsize),j)
    home_ext <- extend(cell, c(max_r_cells, max_r_cells), value=NA)
    home_dist <- distanceFromPoints(home_ext, home[,c("x","y")])
    homerange_kern <- calc(home_dist, fun = function(x){texmex::dgpd(x/1000,
      sigma=scale, xi=shape, u=0)})
    homerange_kern[homerange_kern < min_prob] <- 0
    if (write_homerange == TRUE) {
      k <- home[,id]
      filename <- paste0(output_dir,"/homerange_", k, ".tif")
      writeLines(noquote(paste0("Writing: ", filename)))
      writeRaster(homerange_kern, filename=filename, overwrite=TRUE)
    }
    homerange_kernels[[j]] <- homerange_kern
    names(homerange_kernels[[j]]) <- df_home[j, name]
  }
  return(homerange_kernels)
}

#' Create home range kernels based on Pareto and Gamma functions
#'
#' Creates RasterLayers of homerange kernels based on home centroids
#'
#' @param df_all dataframe of all homerange centroids
#' @param df_home dataframe of homerange centroids to calculate homerange
#'   kernels (these must be a subset of the df_all dataframe), Default is to
#'   use 'df_all' dataframe
#' @param base base Raster that sets the projection, extent, and dimensions of
#'   the study area
#' @param max_r maximum radius to calculate the homerange raster from each
#'   df_home centroid
#' @param pareto_shape shape parameter of the Logistic function that governs
#'   the home kernel
#' @param pareto_scale scale parameter of the Logistic function that governs
#'   the home kernel
#' @param gamma_shape shape parameter of the Gamma function that governs the
#'   conspecific avoidance kernel
#' @param gamma_scale scale parameter of the Gamma function that governs the
#'   conspecific avoidance kernel
#' @param id column name of df_home that identifies the homerange. Default is
#'   NULL, which sets the names to df_home row number.
#' @param output_dir directory for output files (distance, homerange)
#' @param write_distance logical, write distance raster to file. Default is
#'   FALSE.
#' @param write_homerange logical, write home range raster to file. Default is
#'   FALSE
#'
#' @return A list containing homerange kernel Rasters for all the df_home
#'   centroids
#' @export
#'
CreateHomeRangeKernelsParetoGamma <- function(df_all,
                                              df_home = df_all,
                                              base = base,
                                              max_r,
                                              pareto_shape,
                                              pareto_scale,
                                              gamma_shape,
                                              gamma_scale,
                                              id = "nest_id",
                                              output_dir = getwd(),
                                              write_distance = TRUE,
                                              write_homerange = TRUE){
  cellsize <- raster::res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  df_sp <- sp::SpatialPointsDataFrame(df_all[,c("x","y")], df_all,
    proj4string=crs(base))
  base <- raster::crop(base, df_sp, snap='out')
  base <- raster::extend(base, c(max_r_cells/2, max_r_cells/2), value=11)
  writeLines(noquote(paste("Calculating global distance")))
  if (write_distance == TRUE) {
    ifelse(exists("i"), i <- i , i <- 1)
    filename <- paste0(output_dir,"/global_dist_", sprintf("%03d", i), ".tif")
    global_dist <- raster::distanceFromPoints(base, df_sp, filename=filename,
        overwrite=TRUE)
    writeLines(noquote(paste("Writing:", filename)))
  } else {
    global_dist <- raster::distanceFromPoints(base, df_sp)
  }
  homerange_ids <- df_home$nest_id
  total <- length(homerange_ids)
  homerange_kernels <- as.list(setNames(rep(NA,length(homerange_ids)),
    homerange_ids), homerange_ids)
  for (j in 1:nrow(df_home)) {
    writeLines(noquote(paste("Calculating homerange", j, "of", total)))
    home <- df_home[j,]
    home_sp <- SpatialPointsDataFrame(home[,c("x","y")], home,
      proj4string=crs(base))
    xy <- CenterXYInCell(home_sp@coords[,"x"], home_sp@coords[,"y"],
      xmin, ymin, cellsize)
    cell_extent <- extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
      (cellsize/2), xy[2]+(cellsize/2))
    cell <- setValues(raster(cell_extent, crs=projection(base), res=cellsize),j)
    home_ext <- extend(cell, c(max_r_cells, max_r_cells), value=NA)
    home_dist <- distanceFromPoints(home_ext, home[,c("x","y")])
    home_kern <- calc(home_dist, fun = function(x){VGAM::dgpd(x/1000, 0,
      scale=pareto_scale, shape=pareto_shape)})
    global_dist_crop <- crop(global_dist, home_dist)
    cent_dist <- overlay(home_dist, global_dist_crop, fun=function (x,y)
      {ifelse(x != y, NA, x)})
    cent_bounds <- boundaries(cent_dist)
    cent_bounds <- subs(cent_bounds, data.frame(from=c(0,1), to=c(NA,1)))
    edge_dist_abs <- distance(cent_bounds)
    edge_dist <- overlay(cent_dist, edge_dist_abs, fun=function(x,y)
      {ifelse(!is.na(x), y*-1, y*1)})
    edge_dist_shift <- calc(edge_dist, function(x) x - cellStats(edge_dist,
      min))
    avoid_kern <- calc(edge_dist_shift, fun = function(x){dgamma(x/1000, scale =
      gamma_scale, shape = gamma_shape)})
    homerange_kern <- overlay(avoid_kern, home_kern, fun=function(x,y){
      p <- y*x})
    if (write_homerange == TRUE) {
      k <- home[,id]
      filename <- paste0(output_dir,"/homerange_", k, "_",sprintf("%03d",
        i), ".tif")
      writeLines(noquote(paste0("Writing: ", filename)))
      writeRaster(homerange_kern, filename=filename, overwrite=TRUE)
    }
    homerange_kernels[[j]] <- homerange_kern
    names(homerange_kernels[[j]]) <- df_home[j,"nest_id"]
  }
  return(homerange_kernels)
}



#' Exports telemetry data into 'Individuals' folders
#'
#' Splits data by ID and exports .csv files to the Individuals folders
#'
#' @param data Dataframe of location data
#' @param output_dir String, folder for 'individual' files. Default is:
#'   "C:/Work/R/Data/BAEA/Telemetry/Individuals/CSVs"
#' @param output_dir_all String, folder for 'all' files. Default is:
#'   "C:/Work/R/Data/BAEA/Telemetry/Individuals/CSVs"
#' @param tz String, time zone. Default is "Etc/GMT+5"
#'
#' @return Exports .csv files.
#' @importFrom magrittr "%>%"
#' @export
#'
ExportTelemetryDataById <- function(data,
                                    output_dir = file.path("C:/Work/R/Data",
                                      "BAEA/Telemetry/Individuals/CSVs"),
                                    output_dir_all = file.path("C:/Work/R/Data",
                                      "BAEA/Telemetry/All/CSVs"),
                                    tz = "Etc/GMT+5"){
  data <- data
  # All data
  filenames <- list.files(output_dir_all, full.names=TRUE)
  csv_file <- filenames[grep("All", filenames)]
  file.remove(csv_file)
  first_date <- format(min(data$datetime), "%Y-%m-%d")
  last_date <- format(max(data$datetime), "%Y-%m-%d")
  new_csv_file <- paste0(output_dir_all, "/All","_", first_date, "_", last_date,
    ".csv")
  writing <- paste( "Writing:", new_csv_file)
  writeLines(noquote(writing))
  write.csv(data, new_csv_file, row.names=FALSE)
  # By individual
  filenames <- list.files(output_dir, full.names=TRUE)
  ids <- unique(data$id)
  for (i in  1:length(ids)){
    id <- ids[i]
    df <- data[data$id == id,]
    classes <- unlist(sapply(data, class))
    drops <- c("datetime2")
    classes <- classes[!names(classes) %in% drops]
    names(classes)[which(names(classes) == "datetime1")] <- "datetime"
    files <- filenames[grep(paste0(id, "\\_"), filenames, perl=TRUE)]
    csv_file <- files[file_ext(files) == "csv"]
    classes["alt"] <- "numeric"
    if (length(csv_file) > 0){
      df2 <- read.csv(csv_file, header=TRUE, stringsAsFactors=FALSE,
        colClasses=classes)
      file.remove(csv_file)
      tz(df2$datetime) <- tz
      df <- df %>%
        rbind(df2, df) %>%
        dplyr::distinct() %>%
        dplyr::arrange(id, datetime)
      first_date <- format(min(df$datetime), "%Y-%m-%d")
      last_date <- format(max(df$datetime), "%Y-%m-%d")
      new_csv_file <- paste0(output_dir, "/", id, "_", first_date, "_",
        last_date, ".csv")
      writing <- paste( "Writing:", new_csv_file)
      writeLines(noquote(writing))
      write.csv(df, new_csv_file, row.names=FALSE)
    } else {
      first_date <- format(min(df$datetime), "%Y-%m-%d")
      last_date <- format(max(df$datetime), "%Y-%m-%d")
      new_csv_file <- paste0(output_dir, "/", id, "_", first_date, "_",
        last_date, ".csv")
      write.csv(df, new_csv_file, row.names=FALSE)
      writing <- paste( "Writing:", new_csv_file)
      writeLines(noquote(writing))
    }
  }
}

#' Exports kml files for each year
#'
#' Splits data by ID and exports .kmls to the Individuals folders
#'
#' @usage ExportTelemetryKMLByYear(data, output_dir, update, tz)
#'
#' @param data Dataframe of location data.
#' @param output_dir Output directory, default is:
#'   "C:/Work/R/Data/BAEA/Telemetry/Individuals/CSVs"
#' @param update Logical, whether or not to update only the current year's
#'   files. Default is TRUE.
#' @param tz String, timezone of data. Default is "Etc/GMT+5".
#'
#' @return Exports .kml files
#' @export
#'
ExportTelemetryKMLByYear <- function(data,
                                     output_dir = file.path("C:/Work/R/Data",
                                       "BAEA/Telemetry/Individuals/KMLs"),
                                     update = TRUE,
                                     tz = "Etc/GMT+5"){
  data <- data
  ids <- unique(data$id)
  for (i in  1:length(ids)){
    id <- ids[i]
    df <- data[data$id == id,]
    if(update == TRUE){
      current_year <- year(now())
      df_year <- df[df$year == current_year,]
      if (nrow(df_year) > 1){
        ExportKMLTelemetryBAEA(df_year, file = paste0(id, "_", current_year,
          ".kml"), output_dir = output_dir)
      }
    } else {
      years <- year(first(df$datetime)):year(last(df$datetime))
      for (j in 1:length(years)){
        year <- years[j]
        df_year <- df[df$year == year,]
        ExportKMLTelemetryBAEA(df_year, file = paste0(id, "_", year, ".kml"),
          output_dir = output_dir)
      }
    }
  }
}

#' Extracts flight data
#'
#' Extracts the data associated with 'flight' for each bird, based on speed
#'   data.
#'
#' @usage ExtractFlightData(df)
#'
#' @param df String, dataframe with locations.
#' @param min_speed Numeric, minimum speed to be classified as flight.
#'
#' @return A dataframe of "flight" data
#' @export
#'
#' @examples
#'
ExtractFlightSpeed <- function(df,
                               min_speed = 2) {
  df <- df
  # get data that only fits the "flight" criteria
  flight <- df[which(df$speed>min_speed),]
  return(flight)
}

#' Extract movement data
#'
#' Extracts detected movements (non-stationary between locations) data for each
#'   bird.
#'
#' @usage ExtractMovements(df, by, min_step_length, max_step_length)
#'
#' @param df Dataframe with location data which includes "step_length" and
#'   "step_time" columns.
#' @param by String, column name of unique identifier to split data, default is
#'   "id".
#' @param min_step_length Numeric, threshold distance for movement
#'   classification. Default is 50.
#' @param max_step_time Numeric, threshold difference in time for movement
#'   classification. Default is 20.
#'
#' @return A dataframe of selected movements
#' @export
#'
ExtractMovements <- function(df,
                             by = "id",
                             min_step_length = 50,
                             max_step_time = 20){
  df <- df
  movements <- plyr::ddply(df, by, function(df){
    movements <- df[which(df$step_length >= min_step_length &
    df$step_time <= max_step_time), ]})
  return(movements)
}

#' Extracts roost data
#'
#' Extracts dataframe of sex, arr_diff_min, and dep_diff_min
#'
#' @usage ExtractRoostData(df)
#'
#' @param df Dataframe of location data with "sex", "datetime", "behavior",
#'   "hr_after_sunset", and "hr_before_sunrise" columns.
#'
#' @return a dataframe
#' @details Used in other functions.
#'
ExtractRoostData <- function(df = df){
  arrive <- subset(df, behavior=="arrive")
  arrive$arr_diff_min <- as.integer(difftime(arrive$hr_after_sunset,
    arrive$datetime))
  arrive <- subset(arrive, select=c("sex", "arr_diff_min"))
  arrive$dep_diff_min <- NA
  arrive$roost <- "arrive"
  depart <- subset(df, behavior=="depart")
  depart$dep_diff_min <- as.integer(difftime(depart$datetime,
    depart$hr_before_sunrise))
  depart <- subset(depart, select=c("sex", "dep_diff_min"))
  depart$arr_diff_min <- NA
  depart$roost <- "depart"
  output <- rbind(arrive, depart)
  output$diff_min <- NA
  output$diff_min <- ifelse(!is.na(output$arr_diff_min), output$arr_diff_min,
    output$dep_diff_min)
  return(output)
}

#' Filter location by individual and dates
#'
#' Used to filter full dataset to individual(s) within a specified date range.
#'
#' @param df Dataframe with locations.
#' @param id Column name of unique identifier. Default is "id".
#' @param individual String, individual/s (from id column) to keep, format
#'   should be c(id, id), default is to keep all.
#' @param start Start date filter, default is 1970-01-01.
#' @param end End date filter, default is current date.
#' @param behavior String of specific behaviors to maintain. Default is all.
#'
#' @return A dataframe with subsetted rows
#' @export
#'
#' @details Defaults are specific to my file directories and locations
FilterLocations <- function(df = df,
                            id = "id",
                            individual = NULL,
                            start = NULL,
                            end = NULL,
                            behavior = NULL){
  if (is.null(start) || start == ""){
    start <- "2013-01-01"
  }
  if (is.null(end) || end == ""){
    today <- Sys.Date()
    end <- format(today, format="%Y-%m-%d")
  }
  starts <- paste("Start date: ", start, sep="")
  ends <- paste("End date: ", end, sep="")
  writeLines(noquote(starts))
  writeLines(noquote(ends))
  end <- as.POSIXct(end)
  end <- trunc(end, "days") + 60*60*24
  df <- df[df$datetime >= as.POSIXct(start) & df$datetime <= end,]
  row.names(df) <- NULL
  if(all(is.null(individual)) || individual == ""){
    writeLines(noquote ("All individuals included"))
  } else {
    writeLines(noquote(paste("Filtered to individual(s):", individual, sep="")))
    df <- df[df[,id] %in% individual,]
    row.names(df) <- NULL
  }
  if(is.null(behavior)){
    writeLines(noquote ("All behaviors included"))
  } else {
    writeLines(noquote(paste("Filtered to behavior(s):", behavior, sep="")))
    df <- df[df[,"behavior"] %in% behavior,]
    row.names(df) <- NULL
  }
  return(df)
}

#' Compares two times and returns greater or lesser value
#'
#' Examines two time columns, returns the greater or lesser value and inserts
#'   that value in a result column
#'
#' @usage IfElseTimedateCompare(df, col1, col2, result, default_tz, tz)
#'
#' @param df Dataframe of location data.
#' @param col1 Column name, col1 is compared to col2
#' @param sign String, either ">" or "<" for returning greater than or less than
#'   values, respectively. Default is ">".
#' @param col2 Column name, col2 is compared to col1.
#' @param result String, name of column to insert results. Default is "result".
#' @param default_tz String, timezone that timedate reverts to, based on OS
#'   time.
#' @param tz String, timezone of original data (may be different from
#'   default_tz).
#'
#' @return Original dataframe with a new "result" column
#'
#' @details  Used in RoostArrivalDeparture(). Ensures that result is in the
#'   correct timezone.
#'
IfElseTimedateCompare <- function(df = df,
                                  col1 = "col1",
                                  sign = ">",
                                  col2 = "col2",
                                  result = "result",
                                  default_tz = default_tz,
                                  tz = tz) {
  safe.ifelse <- function(cond, yes, no) {
    structure(ifelse(cond, yes, no), class = class(yes))
  }
  df$col1<-df[,col1]
  df$col2<-df[,col2]
  df$result <- NA
  if (sign == ">"){
    df$result <- safe.ifelse(df$col1 >= df$col2, df$col1, df$col2)
  }
  if (sign == "<"){
    df$result <- safe.ifelse(df$col1 <= df$col2, df$col1, df$col2)
  }
  df$result <- force_tz(df$result, tzone = default_tz)
  df$result <- with_tz(df$result, tzone=tz)
  df$col1 <- NULL
  df$col2 <- NULL
  df_length <- length(df)
  colnames(df)[df_length] <- result
  return(df)
}

#' Finds time value that is not NA
#'
#' Examines two time columns, if one is NA, the other time is returned.
#'
#' @usage IfElseTimedateNA(df, col1, col2, result, default_tz, tz)
#'
#' @param df Dataframe of location data.
#' @param col1 Column name, column checked to if is NA. If not NA, then col1 is
#'   returned.
#' @param col2 Column name, if col1 is NA, then col2 is returned.
#' @param result String, name of column to insert results. Default is "result".
#' @param default_tz String, timezone that timedate reverts to, based on OS
#'   time.
#' @param tz string, timezone of original data (may be different from
#'   default_tz).
#'
#' @return A dataframe with a new 'result' column.
#' @export
#'
#' @details Ensure that result is in the correct timezone. Used in
#'   RoostArrivalDeparture().
#'
IfElseTimedateNA <- function(df = df,
                             col1 = "col1",
                             col2 = "col2",
                             result = "result",
                             default_tz = default_tz,
                             tz = tz) {
  safe.ifelse <- function(cond, yes, no) {
    structure(ifelse(cond, yes, no), class = class(yes))
  }
  df$col1<-df[,col1]
  df$col2<-df[,col2]
  df$result <- NA
  df$result <- safe.ifelse(is.na(df$col1), df$col2, df$col1)
  df$result <- force_tz(df$result, tzone = default_tz)
  df$result <- with_tz(df$result, tzone=tz)
  df$col1 <- NULL
  df$col2 <- NULL
  df_length <- length(df)
  colnames(df)[df_length] <- result
  return(df)
}

#' Imports 'baea' data
#'
#' Imports baea.csv and merges with existing data
#'
#' @usage ImportBAEA(existing, import, tz)
#'
#' @param existing Dataframe, exisiting file to merge with baea. Can be NULL.
#'   Default is "deployed".
#' @param import Logical, whether or not to import baea.csv file. Default is
#'   TRUE.
#' @param tz String, timezone. Default is "Etc/GMT+5".
#'
#' @return Merged baea dataframe is returned. Also, writes new "baea.csv" file.
#' @export
#'
#' @details Directory defaults are specific to my computer
#'
ImportBAEA <- function(existing = deployed,
                       import = TRUE,
                       tz = "Etc/GMT+5") {
  if (import == TRUE) {
    baea <- read.csv(file="C:/Work/R/Data/BAEA/BAEA.csv", header = TRUE,
      stringsAsFactors=FALSE)
    date_cols <- c("date","on_hand","deployed", "end_data",  "failed",
      "removed", "recovered")
    for (i in date_cols) {
      baea[,i] <- as.Date(baea[,i], "%Y-%m-%d", tz=tz)  # they are chars
    }
    datetime_cols <- c("datetime","sunset", "sunrise",  "solarnoon",
      "hr_before_sunrise", "hr_after_sunset")
    for (i in datetime_cols) {
      baea[,i] <- as.POSIXct(baea[,i], tz=tz, usetz=FALSE)
    }
    if (!is.null(existing)) {
      existing <- subset(existing, date > (as.Date(max(baea$date)) - days(3)))
      baea <- subset(baea, date <= (as.Date(max(baea$date)) - days(3)))
        max(baea$date)
        min(existing$date)
    # 3 days are removed from baea to ensure that AddSegmentTimeLength, etc. was
    # done on a full dataset. The baea and existing datasets should not overlap.
      if(!("sunrise" %in% colnames(existing))) {
        existing <- AddSolarTimes(existing)
      }
      existing <- AddStepLengthAndAngles(existing)
      existing <- AddStepTime(existing)
      existing <- AddTimeStepProportion(existing)
      existing <- AddFirstLastDistance(existing)
      baea_full <- rbind(baea, existing)
      baea_full <- unique(baea_full)
      baea_full <- baea_full[with(baea_full,order(id,datetime)),]
      row.names(baea_full) <- NULL
      date <- Sys.Date()
      outfile <- paste("C:/Work/R/Data/BAEA/Archive/BAEA_", date, ".csv",
        sep ="")
      if (!file.exists(outfile)) {
        writeLines(noquote(paste("Merging existing and import")))
        write.csv(baea_full, file=outfile, row.names=FALSE)
        writeLines(noquote(c("Writing: ", outfile, sep = "")))
      }
      write.csv(baea_full, file="C:/Work/R/Data/BAEA/BAEA.csv", row.names=FALSE)
        # rewrites import file
      writeLines(noquote("Writing: \"C:/Work/R/Data/BAEA/BAEA.csv\""))
      baea <- baea_full
    }
  }
  if (import == FALSE) {
  writeLines(noquote(paste("BAEA.csv was NOT imported")))
    if (!is.null(existing)) {
      if(!("sunrise" %in% colnames(existing))) {
        existing <- AddSolarTimes(existing)
      }
      existing <- AddStepLengthAndAngles(existing, by = "id")
      existing <- AddStepTime(existing, by = "id")
      existing <- AddTimeStepProportion(existing)
      existing <- AddFirstLastDistance(existing)
      writeLines(noquote(paste("Coverted existing to baea", sep="")))
      baea <- existing
      } else {
        writeLines(noquote("Nothing imported or converted"))
      }
  }
    return(baea)
}

#' Plots daily behavior proportions as bars
#'
#' Plots proportion of behavioral states during a day as a bar graph, with
#'   adjustable number of breaks.
#'
#' @usage PlotBehaviorProportionBar(df, breaks)
#'
#' @param df Dataframe with "sex", "behavior", and "time_proportion" columns.
#' @param breaks Numeric, number of breaks in daily period.
#'
#' @return Facetted plot of behavior proportion over daily period.
#' @export
#'
#' @details Behavioral colors come from:
#'   "C:/Work/R/Data/BAEA/Models/Behavior_Colors.csv".
#'
PlotBehaviorProportionBar<-function(df = df,
                                    breaks = 20){
  behavior_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/Models/Behavior_Colors.csv", metadata_id="behavior")
  df$behavior <- factor(df$behavior)
  CutProportion <- function(data,breaks=breaks) {
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    k <- cut(data, breaks=brk)
  }
  CutProportionMid <- function(data,breaks=breaks) {
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    mid <- b[1:breaks*2]
    k <- cut(data, breaks=brk)
    mid[k]
  }
  df$bins <- CutProportion(df$time_proportion, breaks)
  df$bins_mid <- factor(CutProportionMid(df$time_proportion,breaks))
  melted <- melt(ddply(df, .(sex, bins_mid),
    function(x){prop.table(table(x$behavior))}))
  names(melted)[names(melted) == 'variable'] <- 'behavior'
  melted$bins_mid <- as.numeric(as.character(melted$bins_mid))
  ggplot(melted, aes(x = bins_mid, y=value, ymax=1, fill= behavior)) +
    facet_grid(~ sex) + geom_bar(stat="identity") +
    scale_fill_manual(values = behavior_colors) +
    scale_x_continuous(breaks=seq(0,1,.1)) +
    theme(panel.margin = unit(1, "lines")) +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x="Daily Period",
    y="Behavior Proportion",
    title="Daily Behavior Distributions")
}

#' Plots daily behavior proportions as lines
#'
#' Plots proportion of behavioral states during a day as a line, with
#'   adjustable number of breaks.
#'
#' @usage PlotBehaviorProportionLine(df, breaks)
#'
#' @param df Dataframe with "sex", "behavior", and "time_proportion" columns.
#' @param breaks Numeric, number of breaks in daily period.
#'
#' @return Facetted plot of behavior proportion over daily period.
#' @export
#'
PlotBehaviorProportionLine <- function(df = df,
                                     breaks = 20){
  df$behavior<-factor(df$behavior)
  source('C:/Work/R/Functions/gps.R')
  behavior_colors <- CreateColorsByMetadata(file=
      "C:/Work/R/Data/BAEA/Models/Behavior_Colors.csv", metadata_id="behavior")
  CutProportion <- function(data,breaks=breaks) {
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    k <- cut(data, breaks=brk)
  }
  CutProportionMid <- function(data,breaks=breaks) {
    b <- seq(0, 1, length=2*breaks+1)
    brk <- b[0:breaks*2+1]
    mid <- b[1:breaks*2]
    k <- cut(data, breaks=brk)
    mid[k]
  }
  df$bins <- CutProportion(df$time_proportion, breaks)
  df$bins_mid <- factor(CutProportionMid(df$time_proportion,breaks))
  melted <- melt(ddply(df, .(sex, bins_mid),
    function(x){prop.table(table(x$behavior))}))
  names(melted)[names(melted) == 'variable'] <- 'behavior'
  melted$bins_mid <- as.numeric(as.character(melted$bins_mid))
  ggplot(melted, aes(x = bins_mid, y=value, ymax=1, group= behavior, color=
    behavior)) +  facet_grid(~ sex) + theme(panel.margin=unit(1, "lines")) +
    scale_color_manual(values=behavior_colors)+
    geom_line(stat="identity", size=1.5) +
    scale_x_continuous(breaks=seq(0,1,.1)) +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) + labs(x="Daily Period",
    y="Behavior Proportion", title="Daily Behavior Distributions")
}

#' Plot step lengths histogram
#'
#' Plots a histogram of the step-lengths for each bird
#'
#' @param df Dataframe with flights.
#' @param id String, column name of unique identifier. Default = "id".
#' @param xlim Numeric, x-value limit. Default is NULL.
#' @param bin_width Numeric, bin size. Default is: x-value range/30
#'
#' @return A plot of step-length distances
#' @export
#'
PlotStepLengths <- function(df,
                            id = "id",
                            xlim = NULL,
                            bin_width = NULL){
  df <- df
  sum_move <- SummarizeSE(df, "step_length", "id", na_rm=TRUE)
  id_colors <- CreateColorsByID(output=TRUE)
  grid <- seq(min(df$step_length, na.rm=TRUE), max(df$step_length, na.rm=TRUE),
    length = 100)
  probs = c(0.5, 0.75, 0.95)
  sumstats <- ddply(df, id, function(df){
    quantile(df$step_length, probs = probs, na.rm=TRUE)
  })
  q_probs<-paste("q",probs, sep="")
  colnames(sumstats)[2:(length(probs)+1)] <- q_probs
  if (is.null(xlim)) xlim <- max(df$step_length, na.rm=TRUE)
  g <- ggplot(df, aes(x=step_length, fill=id)) + facet_wrap( ~ id)  +
    scale_fill_manual(values=id_colors) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    xlab("Step length") + ylab("Density") + scale_x_continuous(limits=c(0,xlim))
  if (is.null(bin_width)) bin_width = xlim/30
  g <- g + geom_bar(aes(y = ..density.., fill=id), colour="black",
    binwidth=bin_width)
  build <- ggplot_build(g)
  sumstats$xmax <- max(build$panel$ranges[[1]]$x.range)
  sumstats$ymax <- max(build$panel$ranges[[1]]$y.range)
  g + geom_vline(data=sumstats, aes(xintercept=q0.5), linetype="longdash",
    size=1, colour="black")  +
  geom_vline(data=sumstats, aes(xintercept=q0.75), linetype="dashed",
    size=1, colour="grey20") +
  geom_vline(data=sumstats, aes(xintercept=q0.95), linetype="dashed",
    size=1, colour="grey30") +
  geom_text(data=sumstats, aes(x=q0.5 + (xmax*0.02), y=ymax*.75,
    label=paste("Median:", "\n", signif(q0.5,3), sep="")),
    color="black", hjust=0) +
  geom_text(data=sumstats, aes(x=q0.75+(xmax*0.02), y=ymax*.5,
    label=paste("75%:", "\n", signif(q0.75,3), sep="")),
    color = "grey20", hjust=0) +
  geom_text(data=sumstats, aes(x=q0.95+(xmax*0.02), y=ymax*.25,
    label=paste("95%:", "\n", signif(q0.95,3), sep="")),
    color = "grey30", hjust=0)
}


#' Plot roost behavior ECDF
#'
#' Plots roost empirical distribution funtion and fitted Weibull cumulative
#'   distribution function
#'
#' @usage PlotRoostECDF(df, pars)
#'
#' @param df Dataframe of location data.
#' @param pars List of simulation parameters with male/female, arrive/depart,
#'  shape/scale.
#'
#' @return Facetted plot with roost data and fitted Weibull distributions
#'
#' @import ggplot2
#' @export
#'
#' @details Empirical distribution function extends to 15 min past the end of
#'   last time.
PlotRoostECDF <- function(df = baea,
                          pars = baea_pars) {
  df <- ExtractRoostData(df)
  df2 <- ExtractRoostPars(pars)
  df[df=="m"]<-"male"
  df2[df2=="m"]<-"male"
  df[df=="f"]<-"female"
  df2[df2=="f"]<-"female"
  df2$diff_min<-.85*(max(df$diff_min)+15)
  df2$lab<-paste("shape: ",round(df2$shape,2))
  df2$lab2<-paste("scale: ", round(df2$scale,2))
  vec <- with(df, seq(0, max(diff_min)+15, length=max(diff_min)+16))
  weibull <- ddply(df2, .(sex,roost), function(df) {
    data.frame(diff_min=vec,
    density = pweibull(vec, shape=df$shape, scale=df$scale))
  })
  df2$density<-.25*(max(weibull$density))
  df2$density2<-.15*(max(weibull$density))
  ggplot(df,aes(x=diff_min)) +
  stat_ecdf(aes(colour=sex), size=1) +
  geom_line(aes(x=diff_min, y=density), data=weibull,
            colour="orange", size=1) +
  geom_text(aes(x=diff_min, y=density, label=lab),
    data=df2) +
  geom_text(aes(x=diff_min, y=density2, label=lab2),
    data=df2) +
  facet_grid(roost ~ sex) +
  theme(plot.title=element_text(size=22)) +
  theme(text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(colour="black")) +
  labs(x='Time Difference from Start/End of Sim Day',
   y='Cumulative Probability Density',
   title="Roost Data with Fitted Weibull Distributions",
   colour= "Sex")
}

#' Plots roost behavior probability distribution
#'
#' Plots roost times and fitted Weibull probability distributions
#'
#' @usage PlotRoostHistogram(df, pars)
#'
#' @param df Dataframe of location data
#' @param pars List of simulation parameters with male/female, arrive/depart,
#'   shape/scale
#' @param binwidth Numeric, bin size. Default is 15.
#'
#' @return Facetted plot with roost data and fitted weibull distributions
#' @export
#'
#' @details Weibull probability function extends 15 min past last time.
#'
PlotRoostHistogram <- function(df,
                               pars,
                               binwidth = 15) {
  df <- ExtractRoostData(df)
  df2 <- ExtractRoostPars(pars)
  df[df=="m"] <- "male"
  df2[df2=="m"] <- "male"
  df[df=="f"] <- "female"
  df2[df2=="f"] <- "female"
  df2$diff_min <- .85*max(df$diff_min)+15
  df2$lab <- paste("shape: ", round(df2$shape,2))
  df2$lab2 <- paste("scale: ", round(df2$scale,2))
  grid <- with(df, seq(0, max(diff_min)+15, length=max(diff_min)+16))
  weibull <- ddply(df2, .(sex,roost), function(df) {
    data.frame(diff_min=grid,
    density = dweibull(grid, shape=df$shape, scale=df$scale))
  })
  df2$density <- max(weibull$density)
  df2$density2 <- .85*df2$density
  ggplot(df) +
  geom_histogram(aes(x=diff_min, y=..density..), binwidth=binwidth) +
  geom_line(aes(x=diff_min, y=density), data=weibull,
    colour = "orange", size=1) +
  geom_text(aes(x=diff_min, y=density, label=lab),
    data=df2) +
  geom_text(aes(x=diff_min, y=density2, label=lab2),
    data=df2) +
  facet_grid(roost ~ sex) +
  theme(plot.title=element_text(size=22)) +
  theme(text=element_text(size=20, colour="black")) +
  theme(axis.text=element_text(colour="black")) +
  labs(x='Time Difference from Start/End of Sim Day',
   y='Probability Density',
   title="Histogram of Roost Data with Fitted Weibull Distributions")
}

#' Updates weekly baea data and send emails
#'
#' Compiles and manages csv files and kml files in order to create individual
#'   .csv and .kml files, weekly .kml files, and send an email with weekly .kml
#'   file.
#'
#' @param data Dataframe of deployed data.
#' @param download_recent Logical, whether or not to download recent GPS data,
#'   default is FALSE
#' @param send_email Logical, whether or not to send an email, default is TRUE.
#' @param date Date within week to update, default is now().
#' @param to String, recipients of email, default: c("erynn.call@maine.gov",
#'   "charlie.todd@maine.gov").
#'
#' @return Creates files in "C:/Work/R/Data/BAEA/Telemetry" and
#' @importFrom magrittr "%>%"
#' @export
#'
#' @details Specific to my Bald Eagle Project.
UpdateWeeklyData <- function(data = data,
                             download_recent = FALSE,
                             send_email = TRUE,
                             date = now(),
                             to = c("erynn.call@maine.gov",
                                   "charlie.todd@maine.gov")){
  date <- as.POSIXct(date)
  data <- data
  if (download_recent == TRUE) {
    # Download recent deployed Data
    DownloadCTT(units="deployed", download="recent")
    deployed_recent <- CompileDownloads (units="deployed", compile="recent")
    data <- ImportUnits(units="deployed", existing=deployed_recent,
      import=TRUE)
  }
  # Update local individual files
  ExportTelemetryDataById(data)
  ExportTelemetryKMLByYear(data, update=TRUE)
  # Copy all data files to Google Drive
  input_dir_all = file.path("C:/Work/R/Data/BAEA/Telemetry/All")
  output_dir_all = file.path("C:/Users/Blake/Google Drive/PhD Program",
    "BAEA Project/Telemetry Data/All")
  csv_file_all <- list.files(file.path(input_dir_all, "CSVs"), full.names=TRUE)
  do.call(file.remove, list(list.files(file.path(output_dir_all, "CSVs"),
    full.names=TRUE)))
  file.copy(csv_file_all, file.path(output_dir_all, "CSVs"))
  # Copy individual files to Google Drive
  input_dir = file.path("C:/Work/R/Data/BAEA/Telemetry/Individuals")
  output_dir = file.path("C:/Users/Blake/Google Drive/PhD Program/BAEA Project",
    "Telemetry Data/Individuals")
  csv_files <- list.files(file.path(input_dir, "CSVs"), full.names=TRUE)
  kml_files <- list.files(file.path(input_dir, "KMLs"), full.names=TRUE)
  do.call(file.remove, list(list.files(file.path(output_dir, "CSVs"),
    full.names=TRUE)))
  file.copy(csv_files, file.path(output_dir, "CSVs"))
  file.copy(kml_files, file.path(output_dir, "KMLs"))
  # Filter data to current week
  start_date <- as.character(floor_date(date-period(1, "day"), "week"))
  end_date <- as.character(ceiling_date(date-period(1, "day"), "week"))
  deployed <- FilterLocations(df=deployed_all, id="id", individual="",
    start=start_date, end=end_date)
  # Export KML of Locations
  ExportKMLTelemetryBAEA(df=deployed, file="BAEA Data.kml")
  kml_file <- "C:/Users/Blake/Desktop/BAEA Data.kml"
  start_date2 <- gsub("-", "", start_date)
  end_date2 <- gsub("-", "", end_date)
  kml_output <- file.path(input_dir_all, "KMLs", paste0("BAEA_Data", "_",
    start_date, "_", end_date, ".kml"))
  file.copy(kml_file, kml_output)
  file.copy(kml_output, file.path(output_dir_all, "KMLs"))
  if (send_email == TRUE) {
    library(rJava)
    library(mailR)
    attachments <- kml_output
    mailr_file <- read.csv("C:/Work/R/Data/MailR/mailR.csv", header=FALSE,
      stringsAsFactors=FALSE)
    body <- paste0("Erynn,", "\n\n", "Attached is a KML file of the available ",
      "Bald Eagle GPS location data for ", start_date, " to ", end_date, ".\n",
      "\n", "Additional locations for this period may be added when new data ",
      "becomes available. For the most complete datasets, please refer to the ",
      "compiled individual data located on Google Drive at:",
      "Telemetry Data/Individuals.", "\n",
      "\n", "Thanks,", "\n", "Blake", "\n\n\n", "Blake Massey", "\n",
      "PhD Student", "\n", "University of Massachusetts", "\n",
      "Department of Environmental Conservation", "\n", "160 Holdsworth Way",
      "\n", "Amherst, MA 01003", "\n", "bhmassey@eco.umass.edu","\n",
      "928-254-9221 (cell)")
    subject <- paste("BAEA GPS Data:", start_date, "to", end_date)
    send.mail(from=mailr_file[1,2], to=to, subject=subject, body=body,
      smtp = list(host.name="smtp.gmail.com", port=465, user.name="blakemassey",
      passwd=mailr_file[1,3], ssl=TRUE), attach.files=attachments,
      authenticate=TRUE, send=TRUE)
  }
}

#' Updates weekly .kml files
#'
#' Updates weekly KML files on local and GDrive folders. All weekly intervals
#'   are set to be from Sunday to Saturday, inclusive. All start and end dates
#'   will be adjusted to the beginning or end, respectively, of that week.
#'
#' @param data Dataframe of location data.
#' @param start String, start date.
#' @param end String, end date.
#' @param update_gdrive Logical, whether to update the files on the GDrive.
#'   Default is TRUE
#'
#' @return Creates .kml files in "C:/Work/R/Data/BAEA/Telemetry" and
#'   "C:/Users/Blake/Google Drive/BAEA Project/Telemetry Data" folders.
#' @export
UpdateWeeklyKMLFiles <- function(data,
                                 start,
                                 end,
                                 update_gdrive = TRUE){
  data <- data
  start_date <- floor_date(as.POSIXct(start), "week")
  end_date <- ceiling_date(as.POSIXct(end), "week")
  diff_days <- difftime(end_date, start_date)
  rep_interval <- as.interval(diff_days, start_date)
  step_period <- as.period(1, "week")
  step_period <- step_period
  step_intervals <- list()
  interval_counter <- 1
  current_start <- int_start(rep_interval)
  current_end <- (current_start + step_period)
  stop_point <- int_end(rep_interval)
  while(current_start < (stop_point)) {
    current_end <- (current_start+step_period)
    step_intervals[[interval_counter]] <- interval(current_start,
      current_end)
    interval_counter <- interval_counter + 1
    current_start <- current_start + step_period
  }
  if(int_end(step_intervals[[length(step_intervals)]]) > stop_point){
    step_intervals[[length(step_intervals)]] <- NULL
  }
  kml_folder = file.path("C:/Work/R/Data/BAEA/Telemetry/All/KMLs")
  kml_files <-list.files(kml_folder, full.names=TRUE)
  for (i in 1:length(step_intervals)){
    step_interval <- step_intervals[[i]]
    start_int <- as.character(int_start(step_interval))
    end_int <- as.character(int_end(step_interval))
    file.remove(kml_files[grep(start_int, kml_files)]) # removes old file
    interval_data <- FilterLocations(df=data, id="id", individual="",
      start=start_int, end=end_int)
    if (nrow(interval_data) > 0){
      end_int2 <- as.character(int_end(step_interval)-period(1, "day"))
      kml_file <- paste0("BAEA_Data", "_", start_int, "_", end_int2, ".kml")
      ExportKMLTelemetryBAEA(df=interval_data, output_dir=kml_folder,
        file=kml_file)
    } else {
      end_int2 <- as.character(int_end(step_interval)-period(1, "day"))
      kml_file <- paste0("BAEA_Data", "_", start_int, "_", end_int2,
        "(empty).kml")
      cat("<kml>\n<Document>\n</Document>\n</kml>", file=file.path(kml_folder,
        kml_file), append=TRUE)
    }
    if (update_gdrive == TRUE) {
      kml_gdrive_folder <- file.path("C:/Users/Blake/Google Drive",
        "PhD Program/BAEA Project/Telemetry Data/All/KMLs")
      kml_gdrive_files <- list.files(kml_gdrive_folder, full.names=TRUE)
      file.remove(kml_gdrive_files[grep(start_int, kml_gdrive_files)])
      file.copy(file.path(kml_folder, kml_file), file.path(kml_gdrive_folder,
        kml_file))
    }
  }
}
