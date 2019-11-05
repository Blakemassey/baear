#' Adds cruise behavior
#'
#' Finds cruise behavior that meet given threshold parameters
#'
#' @usage AddCruiseBehavior(df, min_agl, min_speed, threshold_agl)
#' @param df dataframe
#' @param min_agl distance (in meters) above ground
#' @param min_speed minimum speed of bird
#' @param threshold_agl any locations above this threshold of agl ('above ground
#'     level') are labeled "cruise"
#'
#' @return dataframe with behavior column that has "cruise"
#' @export
#'
#' @details should run AddNestBehavior, AddRoostBehavior, and AddCruiseBehavior
#'   prior to this function
AddCruiseBehavior <- function(df = df,
                              min_agl = 200,
                              min_speed = 5,
                              threshold_agl = 250) {
  df_cruise <- df %>%
    dplyr::mutate(bh_cruise = as.character(NA) ) %>%
    dplyr::mutate(bh_cruise = dplyr::if_else(agl >= min_agl &
        speed >= min_speed, "Cruise", bh_cruise)) %>%
    dplyr::mutate(bh_cruise = dplyr::if_else(agl >= threshold_agl, "Cruise",
      bh_cruise)) %>%
    dplyr::mutate(
      bh_cruise_tf = dplyr::if_else(!is.na(bh_cruise), TRUE, FALSE),
      bh_cruise_seq = sequence(rle(bh_cruise_tf)$lengths) * bh_cruise_tf) %>%
    dplyr::select(-bh_cruise_tf)
  return(df_cruise)
}

#' Adds flight behavior
#'
#' Finds location data that meet given threshold parameters for flights, which
#'   include the start point, mid-flight, and final location of a flight.
#'
#' @usage AddFlightBehavior(df, min_speed, min_step_length, max_step_time,
#'     threshold_agl)
#' @param df dataframe
#' @param min_speed minimum speed of bird
#' @param min_step_length distance (in meters) between locations
#' @param max_step_time maximum time between locations
#' @param threshold_agl any locations above this threshold of agl ('above ground
#'     level') are labeled "cruise"
#'
#' @return dataframe with columns for flight_index, flight_step, and
#'  flight_length
#' @export
#'
#' @details adds 'bh_flight' to dataframe
AddFlightBehavior <- function(df,
                              min_speed = 5,
                              min_step_length = 50,
                              max_step_time = 20,
                              threshold_agl = 100){
  df_flight <- df %>%
    plyr::mutate(bh_flight = as.character(NA)) %>%
    plyr::mutate(bh_flight = dplyr::if_else(speed >= min_speed &
        step_time < max_step_time & step_length > min_step_length, "Flight",
        bh_flight)) %>%
    plyr::mutate(bh_flight = dplyr::if_else(agl >= threshold_agl, "Flight",
      bh_flight))
  return(df_flight)
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
      hr_nest_site <- df_home[j,] %>% dplyr::select(id) %>% dplyr::pull()
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
      cent_dist <- overlay(home_dist, global_dist_crop, fun = function(x,y){
        ifelse(x != y, NA, x)})
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
#' @param distance_threshold distance (in meters) from location to nest to assign
#'   a location as "nest" behavior
#'
#' @return dataframe with 'bh_nest' column that has "nest"
#' @export
#'
#' @details need to run AddNestData() prior to running this
#'
AddNestBehavior <- function(df = df,
                            distance_threshold = 50) {
  df <- df
  df <- df %>%
    mutate(bh_nest = ifelse(nest_dist <= distance_threshold, "Nest", NA))
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
  df$bh_perch <- ifelse(df$speed < max_speed, "perch", df$behavior)
  return(df)
}

#' Adds roost behavior to dataframe
#'
#' Finds roost arrivals and departures that meet given threshold parameters
#'
#' @usage AddRoostBehavior(df, at_roost_distance_threshold, depart_timediff_max,
#'   arrive_timediff_max)
#'
#' @param df dataframe
#' @param at_roost_distance_threshold  numeric max distance away from roost
#' @param depart_timediff_max max diff between start of day and depart
#' @param arrive_timediff_max max diff between end of day and arrival
#'
#' @return dataframe with 'bh_roost' column that has arrive, depart, and roost
#' @export
#'
#' @details Adds 'bh_roost' column
AddRoostBehavior <- function(df,
                             at_roost_distance_threshold = 50,
                             depart_timediff_max = 1000,
                             arrive_timediff_max = 1000){
  df <- df
  df_departure <- df %>% filter(datetime < solarnoon)
  df_arrival <- df %>% filter(datetime > solarnoon)
  depart <- df_departure %>%
    group_by(date, id) %>%
    mutate(away_from_roost =
        if_else(dist_first > at_roost_distance_threshold | bh_nest == "Nest",
      1, 0)) %>%
    mutate(depart_roost = away_from_roost == 1 &
      !duplicated(away_from_roost == 1)) %>%
    mutate(lead_depart = lead(depart_roost, 1)) %>%
    mutate(roost = ifelse(lead_depart == 1, "Depart", NA)) %>%
    mutate(pos_depart = which(roost == "Depart")[1]) %>%
    mutate(pos_all = 1) %>%
    mutate(pos_cumsum = cumsum(pos_all)) %>%
    mutate(roost = ifelse(pos_cumsum < pos_depart, "Roost", roost)) %>%
    ungroup() %>%
    dplyr::select(id, datetime, roost) %>%
    mutate(roost_depart = TRUE)
  arrive <- df_arrival %>%
    group_by(date, id) %>%
    arrange(desc(datetime)) %>%
    mutate(away_from_roost =
        if_else(dist_last > at_roost_distance_threshold | bh_nest == "Nest",
      1, 0)) %>%
    mutate(arrive_roost = away_from_roost == 1 &
      !duplicated(away_from_roost == 1)) %>%
    mutate(lead_arrive = lead(arrive_roost, 1)) %>%
    mutate(roost = ifelse(lead_arrive == 1, "Arrive", NA)) %>%
    mutate(pos_arrive = which(roost == "Arrive")[1]) %>%
    mutate(pos_all = 1) %>%
    mutate(pos_cumsum = cumsum(pos_all)) %>%
    mutate(roost = ifelse(pos_cumsum < pos_arrive, "Roost", roost)) %>%
    arrange(id, datetime) %>%
    ungroup() %>%
    dplyr::select(id, datetime, roost) %>%
    mutate(roost_arrive = TRUE)
  df_roost <- full_join(arrive, depart, by = c("id", "datetime")) %>%
    arrange(id, datetime) %>%
    mutate(bh_roost = coalesce(roost.x, roost.y)) %>%
    dplyr::select(id, datetime, bh_roost)
  df_final <- left_join(df, df_roost, by = c("id", "datetime"))
  return(df_final)
}

#' Adds 'season' to dataframe
#'
#' Adds 'season' column to original dataframe. Seasons start at 3-19, 6-20,
#'   9-21, and 12-20 for "Spring", "Summer", "Fall", and "Winter", respectively.
#'
#' @usage AddSeason(df, date_col)
#'
#' @param df dataframe
#' @param by date column. Default is "date".
#'
#' @return dataframe with 'season' column
#' @export
#'
#'
AddSeason <- function(df,
                      date_col = "date"){
  df <- df
  df$date <- df[,date_col]
  numeric_date <- 100*month(df$date) + day(df$date)
  cuts <- base::cut(numeric_date, breaks = c(0, 0319, 0620, 0921, 1220, 1231))
  levels(cuts) <- c("winter", "spring", "summer", "fall", "winter")
  cuts <- as.character(cuts)
  df$season <- cuts
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
      "Data/Assets/behavior_colors.csv", metadata_id="behavior")
    if (by == "id") by_colors <- CreateColorsByMetadata(file=
      "Data/GPS/GPS_Deployments.csv", metadata_id="deploy_location")
    if (by == "sex") by_colors <- CreateColorsByMetadata(file=
      "Data/Assets/behavior_colors.csv", metadata_id="sex")
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

#' Creates conspecific and nest distance raster for each baea in data
#'
#' @param baea dataframe of baea locations
#' @param nest_set dataframe of nests
#' @param base base Raster that sets the projection, extent, and dimensions of
#'   the study area
#' @param output_dir directory for output files (distance, homerange)
#' @param max_r maximum radius to calculate the homerange raster from each
#'   df_home centroid
#' @param write_con_nest_all logical, write con_nest_all raster to file.
#'   Default is TRUE
#'
#' @return baea
#'
#' @importFrom magrittr "%>%"
#' @export
#'
#' @details Creates folders for every year in output directory
#'
CreateConDistRasters <- function(baea,
                                 nest_set,
                                 base = base,
                                 output_dir = "Output/Analysis/Territorial",
                                 max_r = 3000,
                                 write_con_nest_all = TRUE){
  if (!dir.exists(output_dir)) dir.create(output_dir)
  for (k in sort(unique(baea$year))) dir.create(file.path(output_dir, k),
    showWarnings = FALSE)
  cellsize <- raster::res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  for (i in 1:nrow(nest_set)){
    nest_set[i, "x"] <- CenterXYInCell(nest_set[i, "long_utm"],
      nest_set[i, "lat_utm"], xmin, ymin, cellsize)[1]
    nest_set[i, "y"] <- CenterXYInCell(nest_set[i, "long_utm"],
      nest_set[i, "lat_utm"], xmin, ymin, cellsize)[2]
  }
  nest_set_sf <- sf::st_as_sf(x = nest_set, coords = c("x", "y"), crs = 32619)
  for (i in unique(baea$id)) {
    baea_i <- baea %>% dplyr::filter(id == i)
    for (j in sort(unique(baea_i$year))){
      baea_k <- baea %>% dplyr::filter(year == j)
      nest_k_xy <- baea_k %>% dplyr::slice(1) %>% dplyr::select(nest_long_utm,
        nest_lat_utm) %>% as.vector()
      nest_k_id <- baea_k %>% dplyr::slice(1) %>% dplyr::select(nest_site) %>%
        dplyr::pull()
      home_k_x <- CenterXYInCell(nest_k_xy[1], nest_k_xy[2], xmin, ymin,
        cellsize)[[1]] # home nest long
      home_k_y <- CenterXYInCell(nest_k_xy[1], nest_k_xy[2], xmin, ymin,
        cellsize)[[2]] # home nest lat
      home_k_xy <- tibble::tibble(x = home_k_x, y = home_k_y)
      home_k_sf <- sf::st_as_sf(x = home_k_xy, coords = c("x", "y"), crs =32619)
      cell_extent <- raster::extent(home_k_x - (cellsize/2),
        home_k_x + (cellsize/2), home_k_y - (cellsize/2),
        home_k_y + (cellsize/2))
      cell <- raster::setValues(raster(cell_extent, crs = projection(base),
        res = cellsize), j)
      home_ext <- raster::extend(cell, c(max_r_cells, max_r_cells), value = NA)
      summary(home_ext)
      home_dist <- raster::distance(home_ext)
      filter_quo <- paste0("active_", j, " == TRUE")
      nest_set_sf_j <- nest_set_sf %>%
        seplyr::filter_se(filter_quo) %>%
        dplyr::filter(nest_site != nest_k_id) # conspecific nests
      nests_k <- sf::st_contains(st_as_sfc(bb(home_dist)), nest_set_sf_j)
      nest_set_sf_k <- nest_set_sf_j %>% dplyr::slice(unlist(nests_k))
      con_dist <- raster::distanceFromPoints(home_ext,
        st_coordinates(nest_set_sf_k)) # raster of nests
      # Nearest neighbor nest distance at home nest
      home_con_dist <- raster::extract(con_dist, home_k_xy)
      con_dist_home <- raster::calc(con_dist, function(x){home_con_dist - x})
      con_dist_zero <- raster::calc(con_dist_home, function(x){if_else(x >= 0,
        x, 0)})
      con_nest <- raster::overlay(home_dist, con_dist_zero,
        fun = function(x,y){round(x+y)})
      if (write_con_nest_all == TRUE) {
        filename <- file.path(output_dir, j, paste0("ConNest_", i,
          ".tif"))
        raster::writeRaster(con_dist, filename = filename,
          format = "GTiff", overwrite = TRUE)
        writeLines(noquote(paste("Writing:", filename)))
      }
    }
  }
  raster_files <- list.files("Output/Analysis/Territorial",
    pattern = "^ConNest_+.+tif$", full.names = TRUE, recursive = TRUE)
  con_nest <- list()
  for(i in 1:length(raster_files)){con_nest[[i]] <- raster(raster_files[i])}
  con_nest_all <- do.call(merge, con_nest)
  writeRaster(con_nest_all, filename = file.path("Output/Analysis",
    "Territorial", "ConNest_All.tif"), format = "GTiff", overwrite = TRUE)
  return(con_nest_all)
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

#' Create a dataframe of flight path segments
#'
#' Adds flight path segments, which includes locations before and after "flight"
#' behavior. Segments are sequentially numbered for each individual. Used for
#' visualization and analysis of flight paths. Some locations may be represented
#' twice because they are both the start and end of a path.
#'
#' @param df dataframe
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @return dataframe with columns for 'id', 'behavior', and 'path_seg'.
#' @export
#'
#' @details adds 'path_seg' to dataframe
CreateFlightPathSegments <- function(df){
 df_split <- df %>%
    as.data.frame(.) %>%
    group_by(id) %>%
    mutate(datetime_lead = lead(datetime, 1)) %>%
    mutate(step_time2 = as.integer(difftime(datetime_lead, datetime,
      units = "mins"))) %>%
    mutate(step_time_prev = lag(step_time2, 1, default = 0)) %>%
    mutate(split_time = if_else(step_time_prev > 20, 1, 0)) %>%
    mutate(split_time_grp = cumsum(split_time)) %>%
    group_by(split_time_grp) %>%
    mutate(behavior_lead = lead(behavior, 1)) %>%
    mutate(start_flight = if_else(behavior_lead == "Flight" &
        behavior != "Flight", 1, 0)) %>%
    mutate(behavior_lag = lag(behavior, 1)) %>%
    mutate(end_flight = if_else(behavior_lag == "Flight" &
        behavior != "Flight", 1, 0))  %>%
    mutate(bh_flight_tf = if_else(behavior == "Flight", TRUE, FALSE),
      bh_flight_seq = (sequence(rle(bh_flight_tf)$lengths) * bh_flight_tf)) %>%
    ungroup() %>%
    dplyr::select(-c(datetime_lead, step_time2, split_time, split_time_grp,
      step_time_prev, behavior_lead, behavior_lag, bh_flight_tf))
  df_path <- df_split %>%
    group_by(id) %>%
    filter(start_flight != 0 | end_flight != 0 | bh_flight_seq != 0) %>%
    mutate(datetime_lead = lead(datetime, 1)) %>%
    mutate(step_time2 = as.integer(difftime(datetime_lead, datetime,
      units = "mins"))) %>%
    mutate(step_time_prev = lag(step_time2, 1, default = 0)) %>%
    mutate(split_time = if_else(step_time_prev > 20, 1, 0)) %>%
    mutate(split_time_grp = cumsum(split_time)) %>%
    group_by(id, split_time_grp) %>%
    mutate(first_loc = ifelse(row_number()==1, 1, 0)) %>%
    ungroup(.) %>%
    dplyr::select(-c(datetime_lead, step_time2, split_time, step_time_prev))
  df_dup <- df_path %>%
    filter(start_flight == 1 & end_flight == 1) %>%
    mutate(first_loc = 1)
  df_out <- rbind(df_path, df_dup) %>%
    arrange(id, datetime, first) %>%
    group_by(id) %>%
    mutate(path_seg = cumsum(first_loc)) %>%
    ungroup(.) %>%
    dplyr::select(-c(start_flight, end_flight, first_loc))
  return(df_out)
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

#' Filter Locations by Nest Behavior Criteria
#'
#' @usage FilterByNestCriteria(df, min_daily_nest_dist, roll_days,
#'   min_roll_days_nest_dist, seasons)
#'
#' @param df dataframe of locations
#' @param min_daily_nest_dist numeric, bird's daily minimum distance from nest
#'   to meet the criteria
#' @param roll_days numberic, number of days for rolling mean calculation of
#'   nest distance
#' @param min_roll_days_nest_dist numeric, bird's daily minimum distance from
#'   nest to meet the criteria
#' @param seasons seasons to keep (e.g., c("spring","summer")).
#'
#' @return dataframe
#' @export
#'
#' @import dplyr
#'
FilterByNestCriteria <- function(df,
                                 min_daily_nest_dist = 100,
                                 roll_days = NA,
                                 min_roll_days_nest_dist = NA,
                                 seasons = c("spring", "summer")){
  df <- df
  df <- df %>% filter(season %in% seasons)
  df_sum <- df %>%
    group_by(id, date) %>%
    summarize(nest_dist_dy_mean = mean(nest_dist),
      nest_dist_min = min(nest_dist)) %>%
    mutate(nest_dist_met = if_else(nest_dist_min <  min_daily_nest_dist, TRUE,
      FALSE)) %>%
    ungroup()
  if(!is.na(roll_days) && !is.na(min_roll_days_nest_dist)) {
    df_sum <- df %>%
      group_by(id, date) %>%
      mutate(
        nest_dist_run = zoo::rollmean(nest_dist_dy_mean, roll_days, fill=NA,
          na.pad=TRUE, align="left"),
        roll_dist_met = if_else(nest_dist_run <  min_roll_days_nest_dist,
          TRUE, FALSE)) %>%
      ungroup()
  }  else {
    df_sum <- df_sum %>% mutate(roll_dist_met = TRUE)
  }
  df_final <- left_join(df, df_sum, by = c("id", "date")) %>%
      filter(!is.na(nest_dist_met) && !is.na(roll_dist_met)) %>%
      filter(nest_dist_met == TRUE, roll_dist_met == TRUE) #%>%
     # dplyr::select(-c(nest_dist_dy_mean, nest_dist_min, nest_dist_met,
    #    roll_dist_met))
  return(df_final)
}


#' Filter Locations by Roost Criteria
#'
#' @usage FilterByRoostBehavior(df, min_daily_nest_dist, roll_days,
#'   min_roll_days_nest_dist, seasons)
#'
#' @param df dataframe of locations
#' @param overnight_distance_threshold numeric, distance that last PM location
#'   and first AM location must be less than for it to be considered a
#'   confirmed roost location
#' @param number_am_loc numeric, number of GPS locations needed for AM window,
#'   default is 7.
#' @param number_pm_loc numeric, number of GPS locations needed for PM window,
#'   default is 7.
#' @param tz timezone, default is "Etc/GMT+5"
#'
#' @return dataframe
#' @export
#'
#' @import dplyr
#'
FilterByRoostCriteria <- function(df = df,
                                  overnight_distance_threshold = 100,
                                  number_pm_loc = 7,
                                  number_am_loc = 7,
                                  tz = "Etc/GMT+5"){
  df <- df
  df$time_after_start <- as.integer(difftime(df$datetime,
    df$hr_before_sunrise, tz=tz, units = ("mins")))
  df$time_before_end <- as.integer(difftime(df$hr_after_sunset, df$datetime,
    tz=tz, units = ("mins")))
  df$two_hr_after_sunrise <- df$hr_before_sunrise + hours(2)
  df$two_hr_before_sunset <- df$hr_after_sunset - hours(2)
  df$sunrise_window_loc <- df$datetime <= df$two_hr_after_sunrise
  df$sunset_window_loc <- df$datetime >= df$two_hr_before_sunset
  sumstats <- df %>%
    group_by(id, date) %>%
    summarize(
      total_loc = n(),
      am_loc = sum(sunrise_window_loc, na.rm=TRUE),
      pm_loc = sum(sunset_window_loc, na.rm=TRUE)) %>%
    mutate(next_date = lead(date, 1),
      next_day_GPS = (next_date - date),
      next_am_loc = lead(am_loc, 1)) %>%
    ungroup()
    df_sum <- left_join(df, sumstats, by = c("id","date"))
  roost_arrival_confirmed  <- df_sum %>%
    filter(last == "Last" &
        step_length <= overnight_distance_threshold &
        next_day_GPS == 1 &
        pm_loc >= 7 &
        next_am_loc >= 7) %>%
    dplyr::select(id, date)
  roost_departure_confirmed <- roost_arrival_confirmed %>%
    mutate(date = date + 1) # the mornings after "last_roost_confirmed" dates
  roost_arrival_intervals <- inner_join(df, roost_arrival_confirmed,
    by = c("date", "id")) %>% filter(datetime > solarnoon)
    # to "confirm" arrival based on overnight distance
  roost_departure_intervals <- inner_join(df, roost_departure_confirmed,
    by = c("date", "id")) %>% filter(datetime < solarnoon)
  roost_final <- rbind(roost_arrival_intervals, roost_departure_intervals) %>%
    arrange(id, datetime) %>%
    dplyr::select(-c(time_before_end, two_hr_after_sunrise,
      two_hr_before_sunset, sunrise_window_loc, sunset_window_loc))
}


#' Filter locations by individual and dates
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

#' Plots daily behavior proportions as bars
#'
#' Plots proportion of behavioral states during a day as a bar graph, with
#'   adjustable number of breaks.
#'
#' @usage PlotBehaviorProportionBar(df, breaks, title)
#'
#' @param df Dataframe with "sex", "behavior", and "time_proportion" columns.
#' @param breaks Numeric, number of breaks in daily period.
#' @param title Character, main title of plot. Default is "Daily Behavior
#'   Distributions".
#'
#' @return Facetted plot of behavior proportion over daily period.
#' @export
#'
#' @details Behavioral colors come from:
#'   "Data/Assets/behavior_colors.csv".
#'
PlotBehaviorProportionBar <- function(df = df,
                                      breaks = 20,
                                      title = NA){
  if(is.na(title)) title <- "Daily Behavior Distributions"
  behavior_colors <- CreateColorsByMetadata(file=
      "Data/Assets/behavior_colors.csv", metadata_id="behavior")
  df$behavior <- factor(df$behavior)
  CutProportion <- function(data, breaks = breaks) {
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
  Capitalize <- function(string) {
    substr(string, 1, 1) <- toupper(substr(string, 1, 1))
    string
  }
  sex_names <- list('female'="Female", 'male'="Male")
  sex_labeller <- function(variable, value){
    return(sex_names[value])
  }
  df$bins <- CutProportion(df$time_proportion, breaks)
  df$bins_mid <- factor(CutProportionMid(df$time_proportion, breaks))
  melted <- reshape::melt(plyr::ddply(df, plyr::.(sex, bins_mid),
    function(x){prop.table(table(x$behavior))}))
  names(melted)[names(melted) == 'variable'] <- 'behavior'
  melted$bins_mid <- as.numeric(as.character(melted$bins_mid))
  ggplot(melted, aes(x = bins_mid, y=value, ymax=1, fill= behavior)) +
    facet_grid(~ sex, labeller = labeller(sex = Capitalize)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = behavior_colors, name = "Behavior") +
    scale_x_continuous(breaks=seq(0, 1, .1)) +
    theme(panel.spacing = unit(1, "lines")) +
    theme(plot.title=element_text(size=20)) +
    theme(text=element_text(size=18, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x="Daily Period",
      y="Behavior Proportion",
      title=title) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())  +
  theme(panel.background = element_rect(fill = "white", color = "black")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))
}

#' Plots daily behavior proportions as lines
#'
#' Plots proportion of behavioral states during a day as a line, with
#'   adjustable number of breaks.
#'
#' @usage PlotBehaviorProportionLine(df, breaks, title)
#'
#' @param df Dataframe with "sex", "behavior", and "time_proportion" columns.
#' @param breaks Numeric, number of breaks in daily period.
#' @param title Character, main title of plot. Default is "Daily Behavior
#'   Distributions".
#'
#' @return Facetted plot of behavior proportion over daily period.
#' @export
#'
PlotBehaviorProportionLine <- function(df = df,
                                       breaks = 20,
                                       title = NA){
  if(is.na(title)) title <- "Daily Behavior Distributions"
  df$behavior <- factor(df$behavior)
  behavior_colors <- CreateColorsByMetadata(file=
      "Data/Assets/behavior_colors.csv", metadata_id="behavior")
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
  Capitalize <- function(string) {
    substr(string, 1, 1) <- toupper(substr(string, 1, 1))
    string
  }
  df$bins <- CutProportion(df$time_proportion, breaks)
  df$bins_mid <- factor(CutProportionMid(df$time_proportion, breaks))
  melted <- reshape::melt(plyr::ddply(df, plyr::.(sex, bins_mid),
    function(x){prop.table(table(x$behavior))}))
  names(melted)[names(melted) == 'variable'] <- 'behavior'
  melted$bins_mid <- as.numeric(as.character(melted$bins_mid))
  ggplot(melted, aes(x = bins_mid, y=value, ymax=1, group= behavior, color=
    behavior)) + geom_line(stat="identity", size=1.5) +
    facet_grid(~ sex, labeller = labeller(sex = Capitalize)) +
    theme(panel.spacing=unit(1, "lines")) +
    scale_color_manual(values=behavior_colors) +
    scale_x_continuous(breaks=seq(0,1, .1)) +
    theme(plot.title=element_text(size=20)) +
    theme(text=element_text(size=18, colour="black")) +
    theme(axis.text=element_text(colour="black")) + labs(x="Daily Period",
    y="Behavior Proportion", title=title) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())  +
  theme(panel.background = element_rect(fill = "white", color = "black")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"))
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


