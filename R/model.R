#' Submodel for behavioral transisitions in BAEA IBM
#'
#' @param sim
#' @param agent_states
#' @param step_data
#' @param step
#'
#' @return
#' @export
#'
BehaviorSubModelBAEA <- function(sim = sim,
                                 agent_states = agent_states,
                                 step_data = step_data,
                                 step = step){
  verbose <- FALSE
  sex <- agent_states$sex
  beta <- as.matrix(sim$pars$classes[[sex]]$constant$fixed$behavior_betas)
  step_row <- which(step_data$datetime == step)
  current_behavior <- as.numeric(step_data[step_row, "behavior"])
  current_time_prop <- as.numeric(step_data[step_row, "time_proportion"])
  next_time_prop <- as.numeric(step_data[step_row + 1, "time_proportion"])
  julian <- yday(step_data[step_row, "datetime"])
  step_data <- as.data.frame(step_data)
  time_prop_second_to_last <- step_data %>%
    mutate(day = lubridate::date(datetime)) %>%
    filter(day == lubridate::date(step)) %>%
    filter(row_number() == (n() - 1)) %>%
    pull(time_proportion)
  gamma <- diag(5)
  g <- beta[1, ]  #  g = state transition probabilities intercepts
  g <- g +
    beta[2, ] * cos(2*pi * (julian/365)) +
    beta[3, ] * sin(2*pi * (julian/365)) +
    beta[4, ] * cos(2*pi * (current_time_prop)) +
    beta[5, ] * sin(2*pi * (current_time_prop))
  gamma[!gamma] <- exp(g) # Beta values in non-diagonal locations in matrix
  gamma2 <- t(gamma) # probabilities for state transitions are now in rows
  gamma3 <- gamma2/apply(gamma2, 1, sum) # rows sum to 1
  if(current_time_prop <= .5 & current_behavior != 5){
    gamma3[, 5] <- 0
    gamma3 <- gamma3/apply(gamma3, 1, sum)
    gamma3[, 5] <- 0
  }
  # next control is new - it forces agent to leave roost after .3 time prop
  if(current_time_prop >= .3 & current_time_prop <= .5 & current_behavior == 5){
    gamma3[, 5] <- 0
    gamma3 <- gamma3/apply(gamma3, 1, sum)
    gamma3[, 5] <- 0
  }
  if(current_time_prop > .5 & current_behavior == 5){
    gamma3[, 1:4] <- 0
    gamma3[, 5] <- 1
  }
  if(current_behavior == 1){ # prevents 1 -> 5 (Cruise to Roost)
    gamma3[, 5] <- 0
    gamma3 <- gamma3/apply(gamma3, 1, sum)
    gamma3[, 5] <- 0
  }
  if(current_behavior == 5){ # prevents 5 -> 1 (Roost to Cruise)
    gamma3[, 1] <- 0
    gamma3 <- gamma3/apply(gamma3, 1, sum)
    gamma3[, 1] <- 0
  }

  # trans prob. given current behavior
  gamma4 <- gamma3
  next_behavior <- sample(1:5, size = 1, prob = gamma4[current_behavior, ])

  if(step_row != nrow(step_data)) { # prevent the creation of an "extra" step
    if(next_time_prop == time_prop_second_to_last){ # Second to last step
      if (next_behavior != 1){ # For next behavior to be not Cruise
        step_data[step_row + 1, "behavior"] <- next_behavior
      } else { # forces selection of non-Cruise behavior
        gamma4[, 1] <- 0
        gamma4 <- gamma4/apply(gamma4, 1, sum)
        gamma4[, 1] <- 0
        next_behavior <- sample(1:5, size = 1, prob = gamma4[current_behavior,])
        step_data[step_row + 1, "behavior"] <- next_behavior
      }
    }
    if(next_time_prop == 1){ # only Nest or Roost for last location of day
      if(verbose) writeLines("End of day")
      if (next_behavior %in% c(3,5)){
        step_data[step_row + 1, "behavior"] <- next_behavior
      } else { # forces selection of Nest or Roost
        gamma4[, c(1,2,4)] <- 0
        gamma4 <- gamma4/apply(gamma4, 1, sum)
        gamma4[, c(1,2,4)] <- 0
        next_behavior <- sample(1:5, size = 1, prob = gamma4[current_behavior,])
        step_data[step_row + 1, "behavior"] <- next_behavior
      }
    }
    if(current_time_prop == 1){
      step_data[step_row + 1, "behavior"] <- current_behavior
    }
    step_data[step_row + 1, "behavior"] <- next_behavior
  }
  return(step_data)
}

#' CreateTimeStepsBAEA
#'
#' Create the time steps within a step interval based on sunrise and sunset.
#' Specifically designed for BAEA project.
#'
#' @usage CreateTimeStepsBAEA(step_interval)
#'
#' @param step_interval = an Interval object from lubridate that represents the
#' current step interval
#' @param sim = 'sim' list
#'
#' @return a list of intervals
#' @export
#'
CreateTimeStepsBAEA <- function(step_interval = step_interval,
                                agent = agent,
                                sim = sim) {
  time_step_period = sim$pars$global$time_step_period
  options(lubridate.verbose=FALSE)
  start_x <- as.numeric(agent$states$start_x)
  start_y <- as.numeric(agent$states$start_y)
  start_xy_utm<- data.frame(x = agent$states$start_x, y = agent$states$start_y)
  crs_utm <- "+proj=utm +zone=19 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  start_pt <- sp::SpatialPoints(start_xy_utm,
    bbox=NULL, proj4string=sp::CRS("+proj=utm +zone=19 +datum=WGS84"))
  xy_coords <- sp::spTransform(start_pt, sp::CRS("+proj=longlat +datum=WGS84"))
  #sp::coordinates(xy_coords)
  #xy_coords <- matrix(c(-67.7778, 44.8012), nrow = 1)
  interval_start_East <- force_tz(lubridate::int_start(step_interval),
    tzone = "US/Eastern")
  interval_end_East <- force_tz(lubridate::int_end(step_interval)-minutes(1),
    tzone = "US/Eastern") #subtracting one minute to make it the same day
  sunrise <- maptools::sunriset(xy_coords, interval_start_East,
    direction = "sunrise", POSIXct.out = TRUE)[1,2]
  sunset <- maptools::sunriset(xy_coords, interval_end_East,
    direction = "sunset", POSIXct.out = TRUE)[1,2]
  step_interval_start <- round_date(sunrise, "minute") - hours(1)
  step_interval_end <- round_date(sunset, "minute") + hours(1)
  end_int <- lubridate::int_end(lubridate::as.interval(time_step_period,
    step_interval_start))
  steps <- append(step_interval_start, end_int)
  while (end_int < step_interval_end) {
    end_int <- end_int + time_step_period
    steps <- append(steps, end_int)
  }
  if (tail(steps, 1) > step_interval_end) {
    steps[length(steps)] <- step_interval_end
  }
  time_step_list <- list()
  for (i in 1:(length(steps)-1)) {
    interval <- lubridate::as.interval(steps[i], steps[i+1])
    interval@tzone  <- "US/Eastern"
    time_step_list[[length(time_step_list)+1]] <- interval
  }
  return(time_step_list)
}

#' CreateTimeStepsinStepIntervalBAEA
#'
#' Create the time steps within a step interval based on sunrise and sunset in
#'   Bangor, ME. Specifically designed for BAEA project.
#'
#' @usage CreateTimeStepsinStepIntervalBAEA(step_interval, sim)
#'
#' @param step_interval = an Interval object from lubridate that represents the
#' current step interval
#' @param sim = 'sim' list
#'
#' @return a list of intervals
#' @export
#'
CreateTimeStepsInStepIntervalBAEA <- function(step_interval = step_interval,
                                              sim = sim) {
  time_step_period = sim$pars$global$time_step_period
  options(lubridate.verbose=FALSE)
  xy_coords <- matrix(c(-67.7778, 44.8012), nrow = 1)
  interval_start_East <- force_tz(lubridate::int_start(step_interval),
    tzone = "US/Eastern")
  interval_end_East <- force_tz(lubridate::int_end(step_interval)-minutes(1),
    tzone = "US/Eastern") #subtracting one minute to make it the same day
  sunrise <- maptools::sunriset(xy_coords, interval_start_East,
    proj4string = sp::CRS("+proj=longlat +datum=WGS84"), direction ="sunrise",
    POSIXct.out= TRUE)[1,2]
  sunset <- maptools::sunriset(xy_coords, interval_end_East,
    proj4string = sp::CRS("+proj=longlat +datum=WGS84"), direction ="sunset",
    POSIXct.out= TRUE)[1,2]
  step_interval_start <- round_date(sunrise, "minute") - hours(2)
  step_interval_end <- round_date(sunset, "minute") + hours(2)
  end_int <- lubridate::int_end(lubridate::as.interval(time_step_period,
    step_interval_start))
  steps <- tibble(start_time = step_interval_start,
    end_time = end_int - seconds(1)) # sub 1 sec prevents overlap
  while (end_int < step_interval_end) {
    start_int <- end_int
    end_int <- end_int + time_step_period
    steps <- steps %>%
      add_row(start_time = start_int,
               end_time = end_int - seconds(1)) # sub 1 sec prevents overlap
  }
  steps[nrow(steps), "end_time"] <- step_interval_end
  time_step_list <- steps %>%
    mutate(step_interval = lubridate::as.interval(start_time, end_time,
      tzone = "US/Eastern")) %>%
    pull(step_interval)
  return(time_step_list)
}

#' LastDailySubModelBAEA
#'
#' Ensures that the last location in a day has a behavior of either "nest" or
#' "roost".
#'
#' @usage LastDailySubModelBAEA(sim, agent_states, step_data, step)
#'
#' @param sim sim object
#' @param agent_states agent_states object
#' @param step_data step_data object
#' @param step step object
#'
#' @return  step_data
#' @export
#'
LastDailySubModelBAEA <- function(sim = sim,
                                  agent_states = agent_states,
                                  step_data = step_data,
                                  step = step) {
  step_row <- which(step_data$datetime == step)
  current_behavior <- as.numeric(step_data[step_row, "behavior"])
  if (current_behavior %in% c(3,5)){
    step_data[step_row + 1, "behavior"] <- current_behavior
  } else {
    overnight_behavior <- sample(c(3,5), 1)
    step_data[step_row, "behavior"] <- overnight_behavior
    step_data[step_row + 1, "behavior"] <- overnight_behavior
  }
  return(step_data)
}

#' Simplify the sim object so that the sim's object size is smaller
#'
#' @param sim
#'
#' @return
#' @export
#'
SimplifySimSpatialBAEA <- function(sim){
    ssf_source <- sim %>%
      pluck('spatial', 'ssf_layers') %>%
      keep(., is.character)
    sim$spatial <- NULL
    sim[["spatial"]][["ssf_layers"]] <- ssf_source
    return(sim)
}

#' Update agent_state object in BAEA IBM
#'
#' @param agent_states
#' @param sim
#' @param init
#'
#' @return
#' @export
#'
UpdateAgentStatesBAEA <- function(agent_states = NULL,
                              sim = sim,
                              init = FALSE) {
  if (init == TRUE) {
    input <- sim$agents$input
    input <- CreateBirthDate(sim)
    input_columns <- colnames(input)
    na_columns <- c("start_datetime", "died")
    all <- list()
    for (i in 1:nrow(input)) {
      states <- list()
      for (j in input_columns) states <- append(states, input[i, j])
      for (k in 1:length(na_columns)) states <- append(states, NA)
      states <- setNames(states, c(input_columns, na_columns))
      agent <- NamedList(states)
      all <- append(all, NamedList(agent))
    }
    sim$agents <- append(sim$agents, NamedList(all))
    return(sim)
  } else {
    agent_states <- agent_states
    return(agent_states)
  }
}

#' UpdateAgentStepDataBAEA
#'
#' Creates and updates agents step data dataframes with timesteps based on
#' starting locations of each agent and dates for the sim so that the timesteps
#' are similar to the duty cycle of the recorded data: 15-minute time steps
#' starting an hour before sunrise and going until an hour after sunset.
#'
#' @usage UpdateAgentStepDataBAEA(step_data, sim, time_step_dfinit)
#'
#' @param step_data = 'step data' dataframe
#' @param sim = 'sim' object
#' @param init = logical, whether or not this is the initation step
#' @param rep_intervals = 'rep_intervals' object
#'
#' @return an 'step_data' dataframe
#' @export
#'
UpdateAgentStepDataBAEA <- function(step_data = NULL,
                                    sim = sim,
                                    init = FALSE,
                                    rep_intervals = rep_intervals) {
  if (init == TRUE) {
    sim_start <- sim$pars$global$sim_start
    all <- sim$agents$all
    for (i in 1:length(all)) {
      agent <- all[[i]]
      writeLines(paste("Creating initial 'step_data' dataframe for", i, "of",
        length(all)))
      all_time_steps <- as.POSIXct(NA)
      for (j in 1:length(rep_intervals)){
        #j <- 1
        step_intervals <- CreateStepIntervals(rep_intervals[[j]],
          step_period = sim$pars$global$step_period)
        for (k in 1:length(step_intervals)){
          #k <- 1
          time_steps <- CreateTimeStepsBAEA(step_intervals[[k]], agent = agent,
            sim = sim)
          for (m in 1:length(time_steps)){
            all_time_steps <- append(all_time_steps, int_start(time_steps[[m]]))
            if (m == length(time_steps)) all_time_steps <-append(all_time_steps,
              int_end(time_steps[[m]]))
          }
        }
      }   # this is all about getting 'all_time_steps'
      if(is.na(all_time_steps[1])) all_time_steps <- all_time_steps[-1]
      time_steps_df <- data.frame(datetime = all_time_steps) %>%
        mutate(julian = yday(datetime)) %>%
        group_by(julian) %>%
        mutate(day_start = min(datetime),
          day_end = max(datetime),
          day_minutes = as.integer(difftime(day_end,day_start,units="mins")))%>%
        ungroup() %>%
        mutate(time_after_start = as.integer(difftime(datetime, day_start,
          units="mins"))) %>%
        mutate(time_proportion = time_after_start/day_minutes) %>%
        dplyr::select(datetime, julian, time_proportion)
      step_data <- time_steps_df %>%
        mutate(id=agent$states$id,
          behavior = NA,
          x = NA,
          y = NA,
          exp_angle = NA,
          abs_angle = NA) %>%
        dplyr::select(id, datetime, behavior, x, y, exp_angle, abs_angle,
          julian, time_proportion) %>%
        as.data.frame() # when saved as a tibble, broke BehaviorSubModel
      step_data[1, "behavior"] <- 3
      step_data[1, "x"] <- agent$states$start_x
      step_data[1, "y"] <- agent$states$start_y
      agent  <- append(agent, NamedList(step_data))
      all[[i]] <- agent
    }
    sim$agents$all <- all
    return(sim)
  } else {
    step_data <- step_data
    return(step_data)
  }
}

#' UpdateSpatialBAEA
#'
#' Creates and updates population interval dataframe
#'
#' @usage UpdateSpatialBAEA(sim, init)
#'
#' @param sim = sim list, MUST be included in function's parameter argument
#' needed when: when init == TRUE
#' @param init = logical, whether or not this is the initation step
#'
#' @return an 'agents' list
#' @export
#'
UpdateSpatialBAEA <- function(sim = sim,
                              init = TRUE) {
  sim <- sim
  if (init == TRUE) {
    behavior_levels <- c("Cruise", "Flight", "Nest", "Perch", "Roost")
    # male and female same for now
    move_pars <- sim$pars$classes$male$constant$fixed$move_pars %>%
      mutate(
        behavior_num = as.numeric(factor(behavior, levels = behavior_levels)),
        behavior_next_num = as.numeric(factor(behavior_next,
          levels = behavior_levels))) %>%
      mutate(ids = paste0(behavior_num, "_", behavior_next_num))
    move_pars_ids <- move_pars$ids
    move_kernels <- as.list(setNames(rep(NA, nrow(move_pars)),
      move_pars_ids), move_pars_ids)
    for (i in 1:nrow(move_pars)){
      move_pars_i <- move_pars[i, ]
      ignore_von_mises <- ifelse(move_pars_i$behavior[1] %in% c("Cruise",
        "Flight"), FALSE, TRUE)
      kernel_i <- CreateMoveKernelWeibullVonMises(
          max_r = NULL,
          cellsize = 30,
          mu1 = move_pars_i$mvm_mu1[1],
          mu2 = move_pars_i$mvm_mu2[1],
          kappa1 = move_pars_i$mvm_kappa1[1],
          kappa2 = move_pars_i$mvm_kappa2[1],
          mix = move_pars_i$mvm_prop[1],
          shape = move_pars_i$weibull_shape[1],
          scale = move_pars_i$weibull_scale[1],
          ignore_von_mises = ignore_von_mises)
      r <- (30*((nrow(kernel_i)-1)/2))+(30/2)
      kernel_raster <- raster::raster(kernel_i, xmn=-r, xmx=r, ymn=-r, ymx=r)
      move_kernels[[i]] <- kernel_raster
      names(move_kernels[[i]]) <- paste0(move_pars_i$behavior_behavior[1])
    }
    male <- NamedList(move_kernels)
    female <- NamedList(move_kernels)
    classes <- NamedList(male, female)
    sim$spatial[["classes"]] <- classes
    return(sim)
  } else {
#   spatial_timer <- UpdateSpatialTimer(spatial_timer)
#   if (spatial_timer == timer_number) UpdateSpatial(); rm(spatial_timer)
#   spatial <- append()
#   spatial <- sim$spatial
#   sim$spatial <- spatial
    return(sim)
  }
}
