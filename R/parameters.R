#' Adds Normal Distribution Density
#'
#' Returns a vector of normal distribution density values based on the input
#'   values and input data's estimated parameters.
#'
#' @usage AddNormDens(x)
#'
#' @param Numeric vector
#'
#' @return A vector
#' @export
#'
#' @details Designed for use in plyr transform functions
#'
AddNormDens <- function(x) {
  m <- mean(x)
  s <- sd(x)
  dnorm(x, m, s)
}

#' Bin one dimension vector by group
#'
#' Creates a "melted" dataframe of frequencies and probabilities for binned
#'   values
#'
#' @usage Bin1DimByGroup(df, var, group, breaks, style, n, midpoints, melt,
#'   signif, ...)
#'
#' @param df Dataframe with 2 columns.
#' @param var Column name of variable to bin
#' @param group Column name of group
#' @param breaks Numeric, vector of bin ranges
#' @param style String, any "breaks" type compatible with classInterval().
#'   Default is "pretty".
#' @param n Numeric, number of bins, as used by classInterval(). Default is 5.
#' @param midpoints Logical, whether to use midpoints for labels. Otherwise,
#'   ranges are used. Default is TRUE.
#' @param melt Logical, whether not to melt() the data
#' @param signif Integer, the number of digits used in formatting the break
#'   numbers. Default is 5.
#' @param ... Arguments to be passed to the functions called in each style used
#'   by classIntervals()
#'
#' @return Two option: If melt = FALSE, returns the original dataframe with
#'   mid-points or ranges for the variable (with suffixes "_mid" or "_bin"). If
#'   melt = TRUE, a table of frequencies ("freq"), and densitities ("dens") for
#'   each group if "group" parameter is used.
#' @export
#'
Bin1DimByGroup <- function (df,
                            var,
                            group = NULL,
                            breaks = NULL,
                            style = "pretty",
                            n = 5,
                            midpoints = TRUE,
                            melt = FALSE,
                            signif = 5,
                            ... ){
  suppressPackageStartupMessages(library(classInt))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(reshape2))
  vec <- df[,var]
  if (!is.null(breaks)){
    break_pts <- breaks
  } else {
    break_pts <- unlist(classInt::classIntervals(vec, style=style, n=n, ...)[2])
  }
  ints = findInterval(vec, break_pts, all.inside = TRUE)
  if (midpoints == TRUE){
    df[, paste0(var, "_mid")] <- (break_pts[ints] + break_pts[ints + 1]) / 2
  } else {
    df[, paste0(var, "_bin")] <- cut(vec, break_pts)
  }
  if (melt == FALSE) {
    return(df)
  } else {
    if (midpoints == TRUE){
      df[,paste0(var)] <- df[,paste0(var, "_mid")]
    } else {
      df[,paste0(var)] <- df[,paste0(var, "_bin")]
    }
    if (!is.null(group)){
      df_sub <-  df[c(which(colnames(df) == var), which(colnames(df) == group))]
    } else {
      df_sub <-  df[which(colnames(df) == var)]
    }
    df_table <- table(df_sub)
    df_melt <- reshape2::melt(df_table)
    colnames(df_melt)[length(colnames(df_melt))] <- "freq"
    if (!is.null(group)){
      df_melt$group_temp <- df_melt[,group]
      df_melt2 <- as.data.frame(df_melt %>% dplyr::group_by(group_temp) %>%
        mutate("dens" = freq/sum(freq)))
      df_melt3 <- df_melt2 %>% dplyr::select(-group_temp)
      df_melt_final <- df_melt3
    } else {
      df_melt[,"dens"] <- df_melt$freq/sum(df_melt$freq)
      df_melt_final <- df_melt
    }
    return(df_melt_final)
  }
}

#' Bin two-dimension vectors
#'
#' Creates a "melted" dataframe of frequencies and probabilities for binned
#'   values in 2 arrays
#'
#' @usage Bin2Dim(df, breaks1, breaks2, midpoints, include_lowest, signif)
#'
#' @param df Dataframe of 2 columns
#' @param breaks1 Vector of bin ranges or any "breaks" type compatible with
#'   hist(), used on x[,1]. Default = "Sturges"
#' @param breaks2 Vector of bin ranges or any "breaks" type compatible with
#'   hist(), used on x[,2]. Default = "Sturges"
#' @param midpoints Logical, whether or not to use midpoints for labels.
#'   Otherwise, ranges are used. Default = TRUE.
#' @param include_lowest Logical, indicating if an value equal to the lowest
#'   'breaks' value should be included.
#' @param signif Integer, the number of digits used in formatting the break
#'   numbers. Default = 5.
#'
#' @return a dataframe of binned mid-points or ranges, frequencies ("freq"),
#'   and probabilities ("prob")
#' @export
#'
Bin2Dim <- function (df,
                     breaks1 = "Sturges",
                     breaks2 = "Sturges",
                     midpoints = TRUE,
                     include_lowest = TRUE,
                     signif = 5){
  histg1 <- hist(df[,1], breaks = breaks1, plot=FALSE)
  histg2 <- hist(df[,2], breaks = breaks2, plot=FALSE)
  brx <- histg1$breaks
  bry <- histg2$breaks
  freq <- table(cut(df[,1], brx, include.lowest=include_lowest, dig.lab=signif),
    cut(df[,2], bry, include.lowest=include_lowest, dig.lab=signif))
  if (midpoints){
    dimnames(freq)[[1]] <- histg1$mids  # put mid-point as dimnames
    dimnames(freq)[[2]] <- histg2$mids  # put mid-point as dimnames
  }
  melt_df <- reshape2::melt(freq)
  colnames(melt_df) <- c(colnames(df)[1], colnames(df)[2], "freq")
  melt_df[,"prob"] <- melt_df$freq/sum(melt_df$freq)
  return(melt_df)
}

#' Create list index
#'
#' Creates a "indexed" list by changing all the values into sequential numbers
#'
#' @param list List of input data
#'
#' @return List with sequential numbers for all of the values
#' @export
#'
#' @details Works with up to 3 nested levels (lists, sublists, subssublists)
CreateListIndex <- function(list){
  index_count <- 0
  for (i in 1:length(list)) {
    sublist  <- list[[i]]
    if (length(sublist) > 1) {
      for (j in 1:length(sublist)) {
        subsublist <- sublist[[j]]
        if (length(subsublist) > 1) {
          for (k in 1:length(subsublist)) {
            index_count <- index_count + 1
            list[[i]][[j]][k] <- index_count
          }
        } else {
          index_count <- index_count + 1
          list[[i]][[j]] <- index_count
        }
     }
     } else {
       index_count <- index_count + 1
       list[[i]] <- index_count
    }
  }
  return(list)
}

#' Create new parameters list
#'
#' Creates a new list of behavior parameters with either NA or index values.
#'
#' @usage CreateNewPars(index, random)
#'
#' @param index Logical, whether or not the parameter values should be numbers
#'   indicating their relative position. Default is TRUE.
#' @param random Logical, whether or not the parameter values should be filled
#'   in using runif(1,1,2), random=TRUE overrides index=TRUE. Default is FALSE.
#'
#' @return A list of newly generated behavior parameters with NA or index
#'   numbers for each of the values
#' @export
#'
#' @details Used to create a new behavior parameters list from scratch. Uses
#'   ExtractParsIndex().
CreateNewPars <- function(index=TRUE,
                          random=FALSE){
  pars <- list(arrive=list(shape=NA, scale=NA),
               depart=list(shape=NA, scale=NA),
               cruise=list(shape1=NA, shape2=NA),
               forage=list(shape1=NA, shape2=NA),
               flight=list(shape=NA),
               loaf=list(shape1=NA, shape2=NA),
               nest=list(shape1=NA, shape2=NA),
               step=list(location=NA, scale=NA, shape=NA),
               territorial=list(shape1=NA, shape2=NA)
               )
  female<-list(female=pars)
  male<-list(male=pars)
  pars<-c(male,female)
  if (random == TRUE){
    pars <- rapply(pars, f=function(x) ifelse(is.na(x),runif(1,1,2),x),
      how="replace" )
  }
  if (index == TRUE && random == FALSE){
    pars <- CreateListIndex(pars)
  }
  return(pars)
}

#' Density of asymetric generalized von Mises
#'
#' Density of an asymetric generalized von Mises distribution. WORK IN
#'   PROGRESS!
#'
#' @usage DensAsymGenVonMises(x, mu1, mu2, kappa1, kappa2)
#'
#' @param x Vector
#' @param mu1 Numeric, primary direction parameter
#' @param mu2 Numeric, secondary direction parameter
#' @param kappa1 Numeric (non-negative), parameter
#' @param kappa2 Numeric (non-negative), parameter
#'
#' @return NLL of general von Mises distribution
#' @export
#'
DensAsymGenVonMises  <- function(x,
                                 mu1,
                                 mu2,
                                 kappa1,
                                 kappa2){
    d = (mu1-mu2)%%pi
    num <- exp(kappa1*cos(x-mu1) + kappa2*cos(2*(x-mu2)) )
    den <- integrate(function(x){exp(kappa1*cos(x) +
      kappa2*cos(2*(x+d)))},
      0,2*pi)$value
    dens <- num/den
    return(dens)
}

#' Density of asymetric generalized von Mises
#'
#' Density function of a generalized von Mises probability distribtion
#'
#' @usage DensGenVonMises(data, mu1, mu2, kappa1, kappa2)
#'
#' @param x Vector
#' @param mu1 Primary direction parameter
#' @param mu2 Secondary direction parameter
#' @param kappa1 Numeric (non-negative), parameter
#' @param kappa2 Numeric (non-negative), parameter
#'
#' @return NLL of general von Mises distribution
#' @export
#'
#' @details Used in NLL_GenVonMises()
DensGenVonMises <- function(x,
                            mu1,
                            mu2,
                            kappa1,
                            kappa2){
#    if (mu2>pi) mu2 <- mu2-(2*pi)
    d = (mu1-mu2)%%pi
    num <- exp(kappa1*cos(x-mu1) + kappa2*cos(2*(x-mu2)) )
    den <- integrate(function(x){exp(kappa1*cos(x) + kappa2*cos(2*(x+d)))},
      0,2*pi)$value
    dens <- num/den
    return(dens)
}

#' Extract parameters by sex
#'
#' Extracts and returns a dataframe of parameter values by sex
#'
#' @param list List of parameters
#' @param sex String, sex to return: "male", "female", or NULL(both), default
#'   is NULL.
#'
#' @return Dataframe of behavior parameters
#' @export
#'
ExtractParsBySex <- function(list,
                             sex = NULL){
  list <- list
  list2 <- sapply(list, FUN = function(X) unlist(X))
  df <- as.data.frame(list2)
  if (!is.null(sex) && sex == "male") df$female <- NULL
  if (!is.null(sex) && sex == "female") df$male <- NULL
  print(df)
  return(df)
}

#' Extract parameters matrix
#'
#' Extracts a dataframe of parameter values with columns for "sex" and
#'   "behavior"
#'
#' @usage ExtractParsMatrix(pars)
#'
#' @param list List of parameters
#'
#' @return Dataframe of behavior parameters
#' @export
#'
#' @details Useful for comparing parameter values among behaviors
#'
ExtractParsMatrix <- function(list){
  df <- data.frame()
  for (i in 1:length(list)){
    sublist  <- list[[i]]
    df_sub <- rbind.fill(lapply(sublist, as.data.frame))
    behavior <- names(sublist)
    sex <- (rep(names(list)[i], length(behavior)))
    df <- rbind(df,cbind(sex, behavior, df_sub))
  }
  print(df)
  return(df)
}

#' Extracts roost parameters
#'
#' Returns df of roost parameters for male/female, arrive/depart, shape/scale
#'
#' @param pars List of simulation parameters with male/female, arrive/depart,
#'   shape/scale
#'
#' @return Dataframe of roost parameters for male/female, arrive/depart,
#'   shape/scale
#' @export
#'
#' @details Used in PlotRoostECDF
#'
ExtractRoostPars <- function(pars = sim_pars){
  sex <- rep(c("f", "m"),each=2)
  roost <- rep(c("depart", "arrive"),times=2)
  shape <- rep(NA,length=4)
  scale <- rep(NA,length=4)
  output <- data.frame(sex, roost, shape, scale, stringsAsFactors=FALSE)
  for (i in 1:nrow(output)){
  row <- output[i,]
  sex <- row[, "sex"]
  p_sex <-ifelse(sex=="f", "female", "male")
  roost <- row [,"roost"]
  par <- paste(roost, "_pars", sep="")
  output[i, "shape"] <-
    pars[[p_sex]][[par]][['shape']]
  output[i, "scale"] <-
    pars[[p_sex]][[par]][['scale']]
  }
  return(output)
}

#' Fits arrival parameter
#'
#' Fits Weibull parameters to roost arrival data
#'
#' @usage FitArrival(data, shape)
#'
#' @param data Dataframe with a "datetime", "hr_after_sunset" and "behavior"
#'   column
#' @param shape Numeric, starting value of shape, default is 10.
#'
#' @return List of arr_scale, arr_shape
#' @export
#'
FitArrival <- function(data = data,
                       shape = 10) {
  data <- subset(data, behavior=="arrive")
  data$arr_diff_min <- difftime(data$hr_after_sunset, data$datetime)
  data$arr_diff_min <- as.integer(as.numeric(data$arr_diff, units = "mins"))
  arr_mean <- mean(data$arr_diff_min) #calculate mean of arrival differences
  arrive_pars <- bbmle::mle2(NLLWeibull, start=list(shape=shape, mu=arr_mean),
    data=list(data=data$arr_diff_min)) #calculate weibull parameters
  arr_scale <- stats::coef(arrive_pars)[2] / gamma(1 + (1/coef(arrive_pars)[1]))
  names(arr_scale) <- "scale"
  arr_shape <- stats::coef(arrive_pars)[1]
  names(arr_shape) <- "shape"
  z <- c(arr_scale, arr_shape)
  return(z)
}

#' Fits behavior parameter
#'
#' Fit Beta parameters to behavior data
#'
#' @usage FitBehavior(df, shape1, shape2)
#'
#' @param data Dataframe with a "datetime", "hr_after_sunset" and "behavior"
#'   column
#' @param shape1 Starting value of shape1, default is 1.
#' @param shape2 Starting value of shape2, default is 1.
#'
#' @return List of (shape1, shape2)
#' @export
#'
FitBehavior <- function(data,
                        shape1 = 1,
                        shape2 = 1) {
  list_pars <- blank_pars
  breaks = 10
  df <- data
  df$sex <- gsub("m", "male", df$sex)
  df$sex <- gsub("f", "female", df$sex)
  df$behavior <- factor(df$behavior)
  behavior_colors <- CreateColorsByBehavior(output=TRUE)
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
  cruise <- melted %>% filter(behavior == "cruise", sex=="male") %>%
    dplyr::select(bins_mid, value)
  cruise_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=cruise$value))
  list_pars[["male"]][["cruise_pars"]][["shape1"]] <- coef(cruise_mle)[1]
  list_pars[["male"]][["cruise_pars"]][["shape2"]] <- coef(cruise_mle)[2]

  forage <- melted %>%
    filter(behavior == "forage", sex=="male") %>% dplyr::select(bins_mid, value)
  forage_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=forage$value))
  list_pars[["male"]][["forage_pars"]][["shape1"]] <- coef(forage_mle)[1]
  list_pars[["male"]][["forage_pars"]][["shape2"]] <- coef(forage_mle)[2]

  loaf <- melted %>%
    filter(behavior == "loaf", sex=="male") %>% dplyr::select(bins_mid, value)
  loaf_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=loaf$value))
  list_pars[["male"]][["loaf_pars"]][["shape1"]] <- coef(loaf_mle)[1]
  list_pars[["male"]][["loaf_pars"]][["shape2"]] <- coef(loaf_mle)[2]

  nest <- melted %>%
    filter(behavior == "nest", sex=="male") %>% dplyr::select(bins_mid, value)
  nest_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=nest$value))
  list_pars[["male"]][["nest_pars"]][["shape1"]] <- coef(nest_mle)[1]
  list_pars[["male"]][["nest_pars"]][["shape2"]] <- coef(nest_mle)[2]

  roost <- melted %>%
    filter(behavior == "roost", sex=="male") %>% dplyr::select(bins_mid, value)
  roost_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=roost$value), optimizer = "nlminb")
  list_pars[["male"]][["roost_pars"]][["shape1"]] <- coef(roost_mle)[1]
  list_pars[["male"]][["roost_pars"]][["shape2"]] <- coef(roost_mle)[2]

  cruise <- melted %>%
    filter(behavior == "cruise", sex=="female") %>% dplyr::select(bins_mid,
      value)
  cruise_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=cruise$value))
  list_pars[["female"]][["cruise_pars"]][["shape1"]] <- coef(cruise_mle)[1]
  list_pars[["female"]][["cruise_pars"]][["shape2"]] <- coef(cruise_mle)[2]

  forage <- melted %>%
    filter(behavior == "forage", sex=="female") %>% dplyr::select(bins_mid,
      value)
  forage_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=forage$value))
  list_pars[["female"]][["forage_pars"]][["shape1"]] <- coef(forage_mle)[1]
  list_pars[["female"]][["forage_pars"]][["shape2"]] <- coef(forage_mle)[2]

  loaf <- melted %>%
    filter(behavior == "loaf", sex=="female") %>% dplyr::select(bins_mid, value)
  loaf_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=loaf$value))
  list_pars[["female"]][["loaf_pars"]][["shape1"]] <- coef(loaf_mle)[1]
  list_pars[["female"]][["loaf_pars"]][["shape2"]] <- coef(loaf_mle)[2]

  nest <- melted %>%
    filter(behavior == "nest", sex=="female") %>%
    dplyr::select(bins_mid, value)
  nest_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=nest$value)) #calculate weibull parameters
  list_pars[["female"]][["nest_pars"]][["shape1"]] <- coef(nest_mle)[1]
  list_pars[["female"]][["nest_pars"]][["shape2"]] <- coef(nest_mle)[2]

  roost <- melted %>%
    filter(behavior == "roost", sex=="female") %>%
    dplyr::select(bins_mid, value)
  roost_mle <- mle2(BetaNLL, start=list(shape1=shape1, shape2=shape2),
    data=list(data=roost$value), optimizer = "nlminb")
  list_pars[["female"]][["roost_pars"]][["shape1"]] <- coef(roost_mle)[1]
  list_pars[["female"]][["roost_pars"]][["shape2"]] <- coef(roost_mle)[2]

  return(list_pars)
}

#' Fits departure parameters
#'
#' Fit Weibull parameters to roost departure data
#'
#' @usage FitDeparture(data, shape)
#'
#' @param data Dataframe with a "datetime", "hr_before_sunrise", and "behavior"
#'   column
#' @param shape Numeric, starting value of shape. Default is 10.
#'
#' @return  List of (dep_scale, dep_shape)
#' @export
#'
#' @examples
#'
FitDeparture <- function(data=data,
                         shape=10) {
  data <- subset(data, behavior=="depart")
  data$dep_diff_min <- difftime(data$datetime, data$hr_before_sunrise)
  data$dep_diff_min <- as.integer(as.numeric(data$dep_diff, units = "mins"))
  data <- subset(data, select = c("id", "date", "datetime", "behavior",
    "dep_diff_min"))
  dep_mean <- mean(data$dep_diff_min) # calculate mean of departure differences
  depart_pars <- bblme::mle2(NLLWeibull, start=list(shape=shape, mu=dep_mean),
    data=list(data=data$dep_diff_min))  # calculate weibull parameters
  dep_scale <- stats::coef(depart_pars)[2] / gamma(1 + (1/coef(depart_pars)[1]))
    # store departure weibull scale
  names(dep_scale)<-"scale"
  dep_shape <- stats::coef(depart_pars)[1] #store arrive weibull shape
  z <- c(dep_scale, dep_shape)
  return(z)
}

#' Fits Gamma parameter
#'
#' Fits a Gamma distribution to data
#'
#' @param df Dataframe
#' @param var Numeric, variable to fit Gamma distribution
#' @param by Column name to subset data
#' @param scale Numeric, starting value of scale, default is 1
#' @param shape Numeric, starting value of shape, default is 1
#'
#' @return Dataframe of "by" variable and parameters
#' @export
#'
FitGammaParsToData <- function(df,
                               var,
                               by = NULL,
                               scale = 1,
                               shape = 1) {
  df <- df
  ifelse(is.null(by), df$by <- "all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, scale=NA, shape=NA,
    stringsAsFactors=FALSE)
  if (is.null(location)) location <- min(df[,var],  na.rm = TRUE)-1
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    pars_i <- bbmle::mle2(NLLGamma, start=list(scale=scale,
      shape=shape), data=list(data=data[,var]),  method="BFGS")
      #calculates Gamma parameters
    pars[which(pars$by==i),"scale"] = as.numeric(stats::coef(pars_i)[1])
    pars[which(pars$by==i),"shape"] = as.numeric(stats::coef(pars_i)[2])
  }
  names(pars)[names(pars) == 'by'] <- names(pars)[1]
  return(pars)
}

#' Fit generalized Von Mises density to array
#'
#' Fit generalized von Mises probability distribution to an array of x values
#'
#' @usage FitGenVonMisesDensityToArray(var, pars)
#'
#' @param var Variable fitted with generalized von Mises distribution
#' @param pars Parameters of generalized von Mises distribution, created by
#'   FitGenVonMisesParsToData()
#'
#' @return Probability densities of a generalized von Mises distribution along
#'   an array of values
#' @export
#'
#' @details The x-axis is set to have 359 values, starting with 1
FitGenVonMisesDensityToArray <- function(var,
                                         pars) {
  x <- seq(0, 2*pi, by=((2*pi)/720))
  dens <- data.frame(by=rep(pars[,1], each=length(x)), x=rep(x,
    time=length(pars[,1])), y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
    dens[which(dens$by==i),"y"] <- DensGenVonMises(x, mu1=pars[which(pars[,1]
      ==i), "mu1"], mu2=pars[which(pars[,1] == i), "mu2"], kappa1=pars[which(
      pars[,1] == i), "kappa1"], kappa2=pars[which(pars[,1] == i), "kappa2"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

#' Fit generalized von Mises parameter to data
#'
#' Fits a generalized von Mises distribution to data
#'
#' @usage FitGenVonMisesParsToData(df, by, var, mu1, mu2, kappa1, kappa2)
#'
#' @param df Dataframe
#' @param var Variable to fit generalized von Mises distribution
#' @param by Column to subset data
#' @param mu1 Primary direction parameter
#' @param mu2 Secondary direction parameter
#' @param kappa1 Numeric (non-negative) parameter of the distribution
#' @param kappa2 Numeric (non-negative) parameter of the distribution
#'
#' @return Dataframe of "by" variable and parameters
#' @export
#'
FitGenVonMisesParsToData <- function(df,
                                     var,
                                     by = NULL,
                                     mu1 = 3*pi/2,
                                     mu2 = pi/2,
                                     kappa1 = .5,
                                     kappa2 = .5) {
  df <- df
  ifelse(is.null(by), df$by<-"all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, mu1=NA, mu2=NA, kappa1=NA, kappa2=NA,
    stringsAsFactors=FALSE)
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    data <- data[!is.na(data)]
    if (length(data) >  1000) data <- sample(data, 800)
    suppressWarnings(pars_i <- mle2(NLLGenVonMises, start=list(mu1=mu1,
      mu2=mu2, kappa1=kappa1, kappa2=kappa2), data=list(data=data)))
    pars[which(pars$by==i), "mu1"] = coef(pars_i)[1]
    pars[which(pars$by==i), "mu2"] = coef(pars_i)[2]
    pars[which(pars$by==i), "kappa1"] = coef(pars_i)[3]
    pars[which(pars$by==i), "kappa2"] = coef(pars_i)[4]
  }
  names(pars)[names(pars) == 'by'] <- names(pars)[1]
  return(pars)
}

#' Fits a Normal density to an array
#'
#' Fit Normal Probability density function to an array of x values
#'
#' @param df Dataframe with data
#' @param var Column name of variable fitted with normal distribution
#' @param pars Parameters of normal distribution, created by
#'   FitNormalParsToData()
#' @param xlim Numeric, sets x-axis limit max, default is: max(df[,var])
#'
#' @return Probability densities of a normal distribution along an array of
#'   values
#' @export
#'
#' @details The x-axis is set to have 500 values, starting with: min(df[,
#'   var])+1
FitNormalDensityToArray <- function(df,
                                    var,
                                    pars,
                                    xlim = NULL) {
  start <- min(df[, var]) + 1
  if (is.null(xlim)){
    end <- max(df[, var])
    x <- seq(start, end, length=500)
  } else {
    x <- seq(start, xlim, length=500)
  }
  dens <- data.frame(by=rep(pars[,1], each=length(x)), x=rep(x,
    time=length(pars[,1])), y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
  dens[which(dens$by == i), "y"] <- stats::dnorm(x, mean=
    pars[which(pars[,1] == i), "mean"], sd=pars[which(pars[,1] == i),"sd"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

#' Fits normal parameters to data
#'
#' Fits a normal distribution to data
#'
#' @usage FitNormalParsToData(df, var, by)
#'
#' @param df Dataframe
#' @param var Column name of variable to fit Wrapped Normal distribution
#' @param by Column name used to subset data
#'
#' @return Dataframe of "by" variable and parameters
#' @export
#'
FitNormalParsToData <- function(df,
                                var,
                                by = NULL) {
  df <- df
  ifelse(is.null(by), df$by<-"all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, mean=NA, sd=NA, stringsAsFactors=FALSE)
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    data <- data[!is.null(data)]
    pars_i <- fitdistr(data[ ,var], "normal")
    pars[which(pars$by==i), "mean"] = as.numeric(pars_i$estimate[["mean"]])
    pars[which(pars$by==i), "sd"] = as.numeric(pars_i$estimate[["sd"]])
  }
  names(pars)[names(pars) == 'by'] <- names(pars)[1]
  return(pars)
}

#' Fit Pareto density to array
#'
#' Fit Pareto probability density function to an array of x values
#'
#' @usage FitParetoDensityToArray(df, var, pars, xlim)
#'
#' @param df Dataframe with data
#' @param var Column name of variable fitted with Pareto distribution
#' @param pars Parameters of Pareto distribution, created by
#'   FitParetoParsToData()
#' @param xlim Numeric, x-axis limit max, default is: max(df[,var])
#'
#' @return Probability densities of a Pareto distribution along an array of
#'   values
#' @export
#'
#' @details The x-axis is set to have 500 values, starting with: min(df[,
#'   var])+1
FitParetoDensityToArray <- function(df,
                                    var,
                                    pars,
                                    xlim = NULL) {
  start <- min(df[, var])+1
  if (is.null(xlim)){
    end <- max(df[, var])
    x <- seq(start, end, length=500)
  } else {
    x <- seq(start, xlim, length=500)
  }
  dens <- data.frame(by=rep(pars[,1], each=length(x)), x=rep(x,
    time=length(pars[,1])), y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
  dens[which(dens$by == i), "y"] <- VGAM::dgpd(x, location=pars[which(pars[,1]
    == i), "location"], scale=pars[which(pars[,1]==i),"scale"], shape=
    pars[which(pars[,1] == i),"shape"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

#' Fits Pareto parameters to data
#'
#' Fits a generalized Pareto distribution to data
#'
#' @usage FitParetoParsToData(df, var, by, location, scale, shape)
#'
#' @param df Dataframe of input data
#' @param var Column name of variable to fit Pareto distribution
#' @param by Column name to subset data
#' @param location Numeric, starting value of location, default is:
#'   min(df$var)-1.
#' @param scale Numeric, starting value of scale, default is 1.
#' @param shape Numeric, starting value of shape, default is 0.
#'
#' @return Dataframe of "by" variable and parameters
#' @export
#'
FitParetoParsToData <- function(df,
                                var,
                                by = NULL,
                                location = NULL,
                                scale = 1,
                                shape = 0) {
  df <- df
  ifelse(is.null(by), df$by <- "all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, location=NA, scale=NA, shape=NA,
    stringsAsFactors=FALSE)
  if (is.null(location)) location <- min(df[,var],  na.rm = TRUE)-1
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    if (is.null(location)) location <- min(data, na.rm = TRUE)-1
    if (location < 0) location <- 0
    pars_i <- bblmle::mle2(NLLPareto, start=list(location=location, scale=scale,
      shape=shape), data=list(data=data[,var]),  method="Nelder-Mead",
      skip.hessian=TRUE) #calculates Pareto parameters
    pars[which(pars$by==i),"location"]=as.numeric(coef(pars_i)[1])
    pars[which(pars$by==i),"scale"]=as.numeric(coef(pars_i)[2])
    pars[which(pars$by==i),"shape"]=as.numeric(coef(pars_i)[3])
  }
  names(pars)[names(pars) == 'by'] <- names(pars)[1]
  return(pars)
}

#' Fit wrapped Cauchy density to array
#'
#' Fit wrapped Cauchy probability density function to an array of x values
#'
#' @usage FitWrappedCauchyDensityToArray(df, var, pars)
#'
#' @param var Column name of variable fitted with Wrapped Cauchy distribution
#' @param pars Parameters of Wrapped Cauchy distribution, created by
#'   FitWrappedCauchyToData()
#'
#' @return Probability densities of a Wrapped Cauchy distribution along an
#'   array of values
#' @export
#'
#' @details The x-axis is set to have 359 values, starting with 1
FitWrappedCauchyDensityToArray <- function(var,
                                           pars) {
  x <- seq(0, 2*pi, by=(2*pi/720))
  dens <- data.frame(by=rep(pars[,1], each=length(x)),
                   x=rep(x, time=length(pars[,1])),
                   y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
    dens[which(dens$by==i),"y"] <- CircStats::dwrpcauchy(x,
      mu=pars[which(pars[,1] == i), "mu"], rho=pars[which(pars[,1] == i),"rho"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

#' Fits a wrapped Cauchy distribution
#'
#' Fits a wrapped Cauchy distribution to data
#'
#' @param df Dataframe of input data
#' @param var Column name of variable to fit Wrapped Cauchy distribution
#' @param by Column name with factor to subset data
#' @param mu Numeric, starting value of mu, default is 1
#' @param rho Numeric, starting value of rho, default is 0
#'
#' @return Dataframe of "by" variable and parameters
#' @export
#'
FitWrappedCauchyParsToData <- function(df,
                                       var,
                                       by = NULL,
                                       mu = 1,
                                       rho = 0) {
  df <- df
  ifelse(is.null(by), df$by <- "all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, mu=NA, rho=NA, stringsAsFactors=FALSE)
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    data <- data[!is.na(data)]
    pars_i <- CircStats::wrpcauchy.ml(data, mu=0, rho=.5, acc=1e-015)
    pars[which(pars$by==i), "mu"] = as.numeric(pars_i[1])
    pars[which(pars$by==i), "rho"] = as.numeric(pars_i[2])
  }
  names(pars)[names(pars) == 'by'] <- names(pars)[1]
  return(pars)
}

#' Fit wrapped Normal density to array
#'
#' Fit wrapped Normal probability density function to an array of x-values
#'
#' FitWrappedNormalDensityToArray(var, pars)
#'
#' @param var Column name of variable fitted with Wrapped Normal distribution
#' @param pars Parameters of Wrapped Normal distribution, created by
#'   FitWrappedNormalToData()
#'
#' @return The probability densities of a Wrapped Normal distribution along an
#'   array of values
#' @export
#'
#' @details The x-axis is set to have 359 values, starting with 1
#'
FitWrappedNormalDensityToArray <- function(var,
                                           pars) {
  x <- seq(0, 2*pi, by=(2*pi/720))
  dens <- data.frame(by=rep(pars[,1], each=length(x)),
                   x=rep(x, time=length(pars[,1])),
                   y=NA, stringsAsFactors=FALSE)
  bys <- unique(pars[,1])
  for (i in bys){
  dens[which(dens$by == i),"y"] <-
    dwrappednormal(x, mu=pars[which(pars[,1] == i),"mu"],
      rho=pars[which(pars[,1]==i),"rho"])
  }
  names(dens)[names(dens) == 'by'] <- names(pars)[1]
  return(dens)
}

#' Fit wrapped Normal parameters to data
#'
#' Fits a Wrapped normal distribution to data
#'
#' @param df Dataframe of input data
#' @param var Column name, variable to fit Wrapped Normal distribution
#' @param by Column name to subset data
#' @param mu Numeric, starting value of mu, default is 1
#' @param rho Numeric, starting value of rho, default is 0
#'
#' @return Dataframe of "by" variable and parameters
#' @export
#'
FitWrappedNormalParsToData <- function(df,
                                       var,
                                       by = NULL,
                                       mu = 1,
                                       rho = 0) {
  df <- df
  ifelse(is.null(by), df$by<-"all", df$by <- df[,by])
  vars <- unique(df$by)
  pars <- data.frame(by=vars, mu=NA, rho=NA, stringsAsFactors=FALSE)
  for (i in vars){
    data  <- subset(df, by == i, select = var, rm.na=TRUE)
    data <- data[!is.null(data)]
    suppressWarnings(pars_i <- circular::mle.wrappednormal(data, mu=NULL,
      rho=NULL))
    pars[which(pars$by==i), "mu"] = as.numeric(pars_i$mu)
    pars[which(pars$by==i), "rho"] = as.numeric(pars_i$rho)
  }
  names(pars)[names(pars) == 'by'] <- names(pars)[1]
  return(pars)
}

#' Fits roost parameters
#'
#' Fits roost arrival and departure Weibull functions
#'
#' @usage FitRoost(df, pars)
#'
#' @param df Dataframe of location data
#' @param pars Input parameters to have roost pars appended onto
#'
#' @return List of parameters
#' @export
#'
#' @details Uses ArriveFit() and DepartFit()
FitRoost <- function(df,
                     pars=NULL){
  df_m <- subset(df, sex=="m")
  df_f <- subset(df, sex=="f")
  ifelse(!is.null(pars), pars <- pars, pars <- list())
  pars[['male']][['arrive_pars']] <- FitArrival(df_m)
  pars[['female']][['arrive_pars']] <- FitArrival(df_f)
  pars[['male']][['depart_pars']] <- FitDeparture(df_m)
  pars[['female']][['depart_pars']] <- FitDeparture(df_f)
  return(pars)
}

#' Fit step length data
#'
#' Fits a generalized Pareto distribution to step-length data
#'
#' @usage FitStepLength(data, location, scale, shape)
#'
#' @param df Dataframe with a "step_length" column
#' @param location Numeric, starting value of location, default is:
#'   min(data$step_length)-1
#' @param scale Numeric, starting value of scale, default is 1.
#' @param shape Numeric, starting value of shape, default is 0.
#' @param pars
#'
#' @return List of (step_location, step_scale, step_shape)
#' @export

FitStepLength <- function(df = df,
                          location = NULL,
                          scale = 1,
                          shape = 0,
                          pars = NULL) {
  if (is.null(pars)){
    pars <- CreateNewPars()
  }
  df <- df
  pars <- pars
  df$sex <- gsub("m", "male", df$sex)
  df$sex <- gsub("f", "female", df$sex)
  df_m<- subset(df, sex=="male")
  df_f<- subset(df, sex=="female")
  if (is.null(location)) location <- min(df$step_length,  na.rm = TRUE)-1
  step_pars_f <- mle2(NLLPareto, start=list(location=location, scale=scale,
    shape=shape), data=list(data=df_f$step_length),  method="Nelder-Mead",
    skip.hessian=TRUE) #calculates Pareto parameters
  step_pars_m <- mle2(NLLPareto, start=list(location=location, scale=scale,
    shape=shape), data=list(data=df_m$step_length),  method="Nelder-Mead",
    skip.hessian=TRUE) #calculates Pareto parameters
  pars$female$step = list(location=as.numeric(coef(step_pars_f)[1]),
                           scale=as.numeric(coef(step_pars_f)[2]),
                           shape = as.numeric(coef(step_pars_f)[3]))
  pars$male$step = list(location=as.numeric(coef(step_pars_m)[1]),
                        scale=as.numeric(coef(step_pars_m)[2]),
                        shape = as.numeric(coef(step_pars_m)[3]))
  return(pars)
}


#' Logistic by inflection
#'
#' Returns the value of the logistic function using inflection and scale as
#'   input parameters
#'
#' @usage LogisticByInflection(x, inflection, scale)
#'
#' @param x Numeric, x value
#' @param inflection Numeric, logistic function inflection point parameter
#' @param scale Numeric, logistic function scale parameter
#'
#' @return The y value
#' @export
#'
LogisticByInflection <- function(x,
                                 inflection,
                                 scale) {
  y <- (1/(exp((-(x - inflection)) / scale) + 1))
  return(y)
}

#' Negative Log-likelihood of a Beta distribution
#'
#' Beta Negative Log-Likelihood function
#'
#' @param data Vector of data
#' @param shape1 Numeric, shape1
#' @param shape2 Numeriv, shape2
#'
#' @return NLL value
#' @export
#'
#' @details Used in BehaviorFit()
NLLBeta <- function(data,
                    shape1,
                    shape2) {
  -sum(dbeta(data, shape=shape1, shape2=shape2, log=TRUE))
}

#' Negative Log-likelihood of a Gamma
#'
#' Gamma Negative Log-Likelihood function
#'
#' @usage NLLGamma(data, scale, shape)
#'
#' @param data Vector of data
#' @param scale Numeric, shape parameter
#' @param shape Numeric, scale parameter
#'
#' @return NLL of Gamma
#' @export
#'
NLLGamma <- function(data,
                     scale,
                     shape) {
  -sum(dgamma(data, scale=scale, shape=shape, log=TRUE))
}

#' Negative log-likelihood of general von Mises
#'
#' Negative log-likelihood of general von Mises function
#'
#' Usage NLLGenVonMises(data, mu1, mu2, kappa1, kappa2)
#'
#' @param data Vector of data
#' @param mu1 Numeric, primary direction paramter
#' @param mu2 Numeric, secondary direction parameter
#' @param kappa1 Numeric, non-negative parameter
#' @param kappa2 Numeric, non-negative parameter
#'
#' @return NLL of general von Mises distribution
#' @export
#'
NLLGenVonMises <- function(data,
                           mu1,
                           mu2,
                           kappa1,
                           kappa2){
  -sum(log(DensGenVonMises(x=data, mu1=mu1, mu2=mu2, kappa1=kappa1,
    kappa2=kappa2)))
}

#' Pareto negative Log-Likelihood function
#'
#' Negative Log-Likelihood of the Pareto function
#'
#' @usage NLLPareto(data, location, scale, shape)
#'
#' @param data Vector of data
#' @param location Numeric, location parameter ("left edge" of the probability
#'   density)
#' @param scale Numeric, scale parameter
#' @param shape Numeric, shape parameter
#'
#' @return NLL of Pareto distribution
#' @export
#'
#' @details uses the generalized Pareto distribution from the VGAM package
NLLPareto <- function(data,
                      location,
                      scale,
                      shape) {
  -sum(VGAM::dgpd(data, location=location, scale=scale, shape=shape, log=TRUE))
}

#' Negative log-likelihood of the Weibull function
#'
#' Negative Log-Likelihood of the Weibull function based on shape and mu
#'   parameters.
#'
#' @param data Vector, data
#' @param shape Numeric, shape parameter
#' @param mu Numeric, mu parameter
#'
#' @return NLL of Weibull distribution
#' @export
#'
#' @details used in DepartFit() and ArriveFit()
NLLWeibull <- function(data,
                       shape,
                       mu) {
  scale <- mu/gamma(1+(1/shape))
  -sum(dweibull(data, shape=shape, scale=scale, log=TRUE))
}

#' NonlinearRangeRescale
#'
#' Performs a nonlinear range rescale according the complement of a Gamma
#'   cumulative density function.
#'
#' @usage NonlinearRangeRescale(x, min, max, lowBound, upBound, pars,
#'   movement_kernel)
#'
#' @param x value to rescale.
#' @param min minimum value of x on its original scale. Default is NULL.
#' @param max maximum value of x on its original scale. Default is NULL.
#' @param lowBound the minimum value of the desired scale. Default is 1.
#' @param upBound the maximum value of the desired scale. Default is NULL.
#' @param pars a list or numeric vector containing the scale and shape
#'   parameters of a Gamma distribution in that order.
#' @param movement_kernel a raster or extent object used to calculate the value
#'   for upBound.
#' @param negative a logical indicating whether or not the final rescaled value
#'   should be multiplied by -1.
#'
#' @return rescaled x
#' @export
#'
#' @details Based on code written by Javan Bauder and Kevin McGarigal. This
#'   function currently only uses a Gamma distribution. Arguments whose defaults
#'   are NULL are calculated internally based on the Gamma parameters provided
#'   and the extent of the movement_kernel raster or extent.
NonlinearRangeRescale <- function(x,
                                  min=NULL,
                                  max=NULL,
                                  lowBound=1,
                                  upBound=NULL,
                                  pars,
                                  movement_kernel,
                                  negative=TRUE)
{

  if(is.null(min)){
    max_distance <- qgamma(0.99,pars[[2]],scale=pars[[1]])
    min <- pgamma(max_distance,pars[[2]],
                  scale=pars[[1]],lower.tail=FALSE)
  }


  if(is.null(max)){
    max <- pgamma(((pars[[2]]-1)*pars[[1]]),pars[[2]],
                  scale=pars[[1]],lower.tail=FALSE)
  }


  if(is.null(upBound)){
    upBound <- sqrt((xmin(movement_kernel)-xmax(movement_kernel))^2+
                      (ymin(movement_kernel)-ymax(movement_kernel))^2)
  }

  # Get predicted y (cummulative gamma) for x
  y_pred <- pgamma(x,pars[[2]],scale=pars[[1]],lower.tail=FALSE)
  rescale <- lowBound + (((y_pred-min)/(max-min)) * (upBound-lowBound))
  if(negative==TRUE){
    rescale <- rescale*-1
  }
  return(rescale)
}



#' Plot Beta cumulative distribution function
#'
#' Plots Beta cumulative probability distribution function using shape1 and
#'   shape2 parameters.
#'
#' @usage PlotBetaCDF(shape1, shape2)
#'
#' @param shape1 Numeric, shape1 parameter.
#' @param shape2 Numeric, shape2 parameter.
#' @param color String, line color, default is "darkgreen".
#'
#' @return Plot of BetaCDF
#'
#' @import ggplot2
#' @export
#'
#' @details The x-axis is set to go to 500.
#'
PlotBetaCDF <- function(shape1 = 5,
                        shape2 = 5,
                        color = "darkgreen") {
  x <- seq(0, 1, length=500)
  y <- pbeta(x, shape=shape1, shape2=shape2)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Beta Distribution (shape1 = ", signif(shape1,3), ", shape2 = ",
    signif(shape2,3), ")",  sep="")
  g <- ggplot(df, aes(x, y)) +
    geom_line(color=color, size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g+labs(x='Value', y='Cumulative Probability Density', title=main)
}

#' Plot Beta probability density function
#'
#' Plots Beta probability density function using shape1 and shape2 parameters.
#'
#' @usage PlotBetaCDF(shape1, shape2)
#'
#' @param shape1 Numeric, shape1 parameter.
#' @param shape2 Numeric, shape2 parameter.
#' @param color String, line color, default is "darkgreen".
#'
#' @return Plot of BetaPDF
#'
#' @import ggplot2
#' @export
#'
#' @details The x-axis is set to go to 500.

PlotBetaPDF <- function(shape1 = 5,
                        shape2 = 5,
                        color = "darkgreen") {
  x <- seq(0, 1, length=500)
  y <- dbeta(x, shape=shape1, shape2=shape2)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Beta Distribution (shape1 = ", signif(shape1,3), ", shape2 = ",
    signif(shape2,3), ")",  sep="")
  g <- ggplot(df,aes(x,y)) +
    geom_line(color=color, size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g+labs(x = 'Value', y= 'Probability Density', title=main)
}

#' Plot Cauchy probability density function
#'
#' Plot Cauchy probability density function using location and scale parameters
#'
#' @param location Numeric, location. Default is 5.
#' @param scale Numeric, scale. Default is 1.
#' @param xlim Numeric, sets axis limits of plot. Default is: c(1st, 99th
#'   quantile)
#' @param color String, line color. Default is "darkgreen".
#'
#' @return Plot of Cauchy PDF
#' @export
#'
#' @details The x-axis is set to have 500 values
#'
PlotCauchyPDF <- function(location = 5,
                          scale = 5,
                          xlim = NULL,
                          color = "darkgreen") {
  if (is.null(xlim)){
    x <- seq(stats::qcauchy(.01, location, scale), stats::qcauchy(.99,
      location, scale), length=500)
  } else {
    x <- seq(xlim[1], xlim[2], length=500)
  }
  y <- stats::dcauchy(x, location=location, scale=scale)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Cauchy Distribution (location = ", signif(location,3),
    ", scale = ", signif(scale,3), ")",  sep="")
  g <- ggplot(df,aes(x,y)) +
    geom_line(color=color, size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x = 'Value', y= 'Probability Density', title=main) +
  geom_vline(xintercept = stats::qcauchy(.95, location, scale), size=1,
    colour = "red", linetype = "longdash") +
  geom_vline(xintercept = stats::qcauchy(.05, location, scale), size=1,
    colour = "red", linetype = "longdash")
}

#' Plots Gamma probability distribution function
#'
#' Plots Gamma probability distribution function using shape and scale
#'  parameters.
#'
#' @param shape Numeric, shape parameter. Default is 1.
#' @param scale Numeric, scale parameter. Default is 1.
#' @param col String, line color. Default is "slateblue4".
#' @param max_x Numeric, maximum valur on x scale. Default is 25.
#'
#' @return Plot of proability distribution
#' @export
#'
PlotGammaPDF <- function(shape = 1,
                         scale = 1,
                         col = "slateblue4",
                         max_x = 25){
  x <- seq(0, max_x, length=101)
  y <- dgamma(x, shape=shape, scale=scale)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Gamma Distribution (scale = ", scale, ", shape = ",
    shape, ")", sep="")
  g <- ggplot(df, aes(x, y)) +
    geom_line(colour=col, size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x='Value', y='Probability Density', title=main)
}

#' Plot histogram of data and generalized von Mises distribution
#'
#' Plots a histogram of the data with an option to fit a generalized von Mises
#'   distribution
#'
#' @usage PlotHistAndGenVonMises (df, var, by, pars, bin_width,
#'   fit_gen_von_mises, fit_color, x_lab, title)
#'
#' @param df Dataframe of data
#' @param var Column name, variable to fit von Mises distribution
#' @param by Column name, optional, column name used to subset data, default
#'   is NULL.
#' @param pars Parameter valuess, optional. A set of parameters used in place
#'   of pars generated within the function. Default is NULL.
#' @param bin_width Numeric, bin size, default is: 15 degrees or 2*pi/24 radians
#' @param fit_gen_von_mises Logical, whether or not to fit and show
#'  generalized von Mises distribution. Default is TRUE.
#' @param fit_color String, color used for von Mises fit line, quantiles, and
#'   parameter value text. Default is "black"
#' @param x_lab String, name for x-axis, default is 'var'.
#' @param title String, title of plot, default is NULL.
#'
#' @return Plots of the data with a fitted von Mises distribution.
#' @export
#'
#' @details Automatically adjusts plot for degrees or radians input, but all
#'   parameter estimates are based on radians
PlotHistAndGenVonMises <- function(df,
                                   var,
                                   by = NULL,
                                   pars = NULL,
                                   bin_width = NULL,
                                   fit_gen_von_mises = TRUE,
                                   fit_color = "black",
                                   x_lab = NULL,
                                   title = NULL) {
  ifelse(max(df[,var], na.rm=TRUE) <= 2*pi, radian <- TRUE, radian <- FALSE)
  if (is.null(x_lab)) x_lab <- var
  if (is.null(bin_width)) bin_width = (2*pi)/24
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  if(is.null(by)){
    df$by <- "all"
    by <- "by"
  } else {
     df$by <- df[,by]
  }
  df$var <- df[,var]
  by_colors <- suppressWarnings(CreateColorsByAny(by=by, df=df))
  breaks <- seq(0, (2*pi), by=((2*pi)/12))
  breaks <- breaks[-length(breaks)]
  minor_breaks <- seq(0, 2*pi, by=bin_width)
  limits <- c(0, 2*pi)
  labels <- round(breaks, 2)
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    theme(axis.text.x = element_text(colour="grey20",size=12))+
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(vjust=3)) +
    xlab(x_lab) + ylab("Density") + labs(title = title)
  g <- g + geom_bar(aes(y = ..density.., fill=by), color="black",
    binwidth=bin_width) + coord_polar(start=(1.5*pi), direction = -1) +
    scale_x_continuous(limits=limits, breaks=breaks, minor_breaks=minor_breaks,
      labels = labels) + scale_y_continuous(labels = NULL)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  if (fit_gen_von_mises == TRUE) {
    if (is.null(pars)) {
      pars <- FitGenVonMisesParsToData(df=df, var="var", by="by")
    if(is.null(by)){
        pars$by <- "all"
      } else {
        pars[,length(pars)+1] <- pars$by
        colnames(pars)[length(pars)] <- by
      }
    }
    dens <- suppressWarnings(FitGenVonMisesDensityToArray(var=var, pars=pars))
    if(is.null(by)){
      dens$by <- "all"
    } else {
      dens[,length(dens)+1] <- dens$by
      colnames(dens)[length(dens)] <- by
    }
    g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
    build <- ggplot_build(g)
    pars$xmax <- max(build$panel$ranges[[1]]$theta.range)  # for geom_text
    pars$ymax <- max(build$panel$ranges[[1]]$r.range)  # for geom_text
    g <- g + geom_text(data = pars, aes(y = ymax*1, x = xmax*.365,
      label = paste0("mu1: ", signif(mu1, 3), "\n", "mu2: ", signif(mu2, 3))),
      size=4.5, color = fit_color, hjust=1, vjust=-1.3)
    g <- g + geom_text(data = pars, aes(y = ymax*1, x = xmax*.135,
      label = paste0("kappa1: ", signif(kappa1, 2), "\n", "kappa2: ",
      signif(kappa2, 2))), size=4.5, color = fit_color, hjust=0, vjust=-1.3)
  }
  g
}

#' Plots histogram and Normal distribution
#'
#' Plots a histogram of the data with an option to fit a Normal distribution.
#'
#' @usage PlotHistAndNormal(df, var, by, pars, xlim, bin_width, fit_normal,
#'   fit_color, x_lab, title, hold_axes, labels, lines)
#'
#' @param df Dataframe of data
#' @param var Column name of variable to fit Normal distribution
#' @param by Column name used to subset data, optional. Default is NULL
#' @param pars Parameters, optional. A set of parameters used in place of pars
#'   generated within the function. Default is NULL.
#' @param xlim Numeric, x-value axis limits, if only one value is given, that
#'   is used as the max limit. Default is NULL.
#' @param bin_width Numeric, bin size. Default is: x-value range/30
#' @param fit_normal Logical, whether or not to fit and show Normal
#'   distribution. Default is TRUE.
#' @param fit_color String, color used for Normal fit line, quantiles, and
#'   parameter value text. Default is "orangered".
#' @param x_lab String, label for x-axis. Default is variable name.
#' @param title String, main title name. Default is NULL.
#' @param hold_axes Logical, hold axes even if parameter distribution goes off
#'   of the screen. Default is TRUE.
#' @param labels Logical, show quantile labels. Default is TRUE.
#' @param lines Logical, show quantile lines. Default is TRUE.
#'
#' @return Plots with a Normal distribution fitted for each 'by' factor
#' @export
#'
PlotHistAndNormal <- function(df,
                              var,
                              by = NULL,
                              pars = NULL,
                              xlim = NULL,
                              bin_width = NULL,
                              fit_normal = TRUE,
                              fit_color = "orangered",
                              x_lab = NULL,
                              title = NULL,
                              hold_axes = TRUE,
                              labels = TRUE,
                              lines = TRUE) {
  if (length(xlim) == 2) xlim <- c(min(xlim), max(xlim))
  if (length(xlim) == 1) xlim <- c(min(df[,var], na.rm=TRUE), xlim)
  if (is.null(xlim)) xlim <- c(min(df[,var], na.rm=TRUE),
    max(df[, var], na.rm=TRUE))
  if (is.null(bin_width)) bin_width = diff(xlim)/30
  if (is.null(x_lab)) x_lab <- var
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])
  df$var <- df[,var]
  by_colors <- CreateColorsByAny(by=by, df=df)
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    xlab(x_lab) + ylab("Density") + labs(title = title) +
    scale_x_continuous(limits=xlim) +
  geom_bar(aes(y = ..density.., fill=by), color="black", binwidth=bin_width)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  ## For "hold_axes"
  if (hold_axes == TRUE) {
    build <- ggplot_build(g)
    g <- g + coord_cartesian(ylim= c(min(build$panel$ranges[[1]]$y.range),
      max(build$panel$ranges[[1]]$y.range)))
  }
  ## Create and plot 'pars' (if fit_normal = TRUE)
  if(fit_normal == TRUE){
    if(is.null(pars)) pars <- FitNormalParsToData(df, var="var", by="by")
    if(is.null(by)){
      pars$by <- "all"
    } else {
      pars[,length(pars)+1] <- pars$by
      colnames(pars)[length(pars)] <- by
    }
    dens <- FitNormalDensityToArray(df=df, var=var, pars=pars)
    if(is.null(by)){
      dens$by <- "all"
    } else {
      dens[,length(dens)+1] <- dens$by
      colnames(dens)[length(dens)] <- by
    }
    g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
  }
  ## Set axes for the rest of the plots
  build <- ggplot_build(g)
  xmax <- max(build$panel$ranges[[1]]$x.range)
  ymax <- max(build$panel$ranges[[1]]$y.range)
  ## Create and plot 'quantiles'
  probs = c(0.05, 0.5, 0.95)
  quantiles <- ddply(df, "by", function(df) quantile(df$va,probs=probs,
    na.rm=TRUE))
  colnames(quantiles)[2:(length(probs)+1)] <- paste("q",probs, sep="")
  quantiles$xmax <- xlim[1]+((xlim[2]-xlim[1])*.975)  # for geom_text
  quantiles$ymax <- ymax  # for geom_text
  if (labels == TRUE) {
  g <- g + geom_text(data = quantiles,
      aes(x = c(q0.05+(abs(xmax)*0.05), q0.5+(abs(xmax)*.05),
          q0.95+(abs(xmax)*.05)),
        y = c(ymax*.8, ymax*.8, ymax*.8),
        label = c(paste("5th:","\n", signif(q0.05, 3)),
        paste("Median:","\n", signif(q0.5, 3)),
        paste("95th:","\n", signif(q0.95, 3)))),
      color = c("black", "gray20", "gray30"), hjust=0, vjust=1)
  }
  if (lines == TRUE) {
  g <- g + geom_vline(data = quantiles, aes(xintercept = c(q0.05, q0.5, q0.95)),
      color = c("black", "gray20", "gray30"), size = 1,
      linetype = c("dashed", "longdash", "dashed"))
  }
  ## Plot 'pars' quantile lines and parameter value text
  if(fit_normal == TRUE) {
    pars$xmax <- xlim[1]+((xlim[2]-xlim[1])*.975)  # for geom_text
    pars$ymax <- ymax   # for geom_text
    if (labels == TRUE) {
    g <- g + geom_text(data = pars, aes(x = xmax, y = ymax*.95,
      label = paste("Mean: ", signif(mean,3),"\n",
      "SD: ", signif(sd,3)), sep=""),
      color = fit_color, hjust=1, vjust=1)
    }
    if (lines == TRUE) {
    g <- g + geom_vline(data = pars, aes(xintercept =
      c(qnorm(.05, mean=mean, sd=sd),
        qnorm(.5, mean=mean, sd=sd),
        qnorm(.95, mean=mean, sd=sd))),
      linetype=c("dashed", "longdash", "dashed"), colour=fit_color, size=1)
    }
    if (labels == TRUE) {
    g <- g + geom_text(data = pars,
      aes(x = c((qnorm(.05, mean=mean,sd=sd) + (abs(xmax)*0.05)),
        (qnorm(.5, mean=mean, sd=sd) + (abs(xmax)*0.05)),
        (qnorm(.95, mean=mean, sd=sd)) + (abs(xmax)*0.05)),
         y = c(ymax*.9, ymax*.9, ymax*.9),
      label = c(paste("5th:","\n",signif(qnorm(.05, mean=mean, sd=sd),3)),
      paste("Median:","\n",signif(qnorm(.5, mean=mean, sd=sd),3)),
      paste("95th:","\n",signif(qnorm(.95, mean=mean, sd=sd),3)))),
      color= fit_color,hjust=0, vjust=1)
    }
  }
  g
}

#' Plots a histogram of the data and a Pareto distribution
#'
#' Plots a histogram of the data with an option to fit a Pareto distribution
#'
#' @usage PlotHistAndPareto(df, var, by, pars, xlim, bin_width, fit_pareto,
#'   fit_color, x_lab, title, hold_axes, emp_labels, emp_lines, labels, lines)
#'
#' @param df Dataframe of data
#' @param var Column name of variable to fit Pareto distribution
#' @param by Column name used to subset data, optional. Default is NULL.
#' @param pars Parameters used in place of pars generated within the function.
#'   Optional. Default is NULL.
#' @param xlim Numeric, x-value limit. Default is NULL.
#' @param bin_width Numeric, bin size. Default is x-value range/30.
#' @param fit_pareto Logical, whether or not to fit and show Pareto
#'   distribution. Default is TRUE.
#' @param fit_color String, color used for Pareto fit line, quantiles, and
#'   parameter value text. Default is "orangered".
#' @param x_lab String, x-axis label. Default is NULL.
#' @param title String, title of plot. Default is NULL.
#' @param hold_axes Logical, hold axes even if parameter distribution goes off
#'   of the screen. Default is TRUE.
#' @param emp_labels Logical, whether or not to show labels for empirical data
#'   75th and  95th quantiles. Default is TRUE.
#' @param emp_lines Logical, whether or not to show lines for empirical data
#'   75th and  95th quantiles. Default is TRUE.
#' @param labels Logical, show quantile labels. Default is TRUE.
#' @param lines Logical, show quantile lines. Default is TRUE.
#'
#' @return Plots with a Pareto distribution fitted for each 'by' factor
#' @export
#'
PlotHistAndPareto <- function(df,
                              var,
                              by = NULL,
                              pars = NULL,
                              xlim = NULL,
                              bin_width = NULL,
                              fit_pareto = TRUE,
                              fit_color = "orangered",
                              x_lab = NULL,
                              title = NULL,
                              hold_axes = TRUE,
                              emp_labels = FALSE,
                              emp_lines = FALSE,
                              labels = TRUE,
                              lines = TRUE) {
  if (is.null(xlim)) xlim <- max(df[, var], na.rm=TRUE)
  if (is.null(bin_width)) bin_width = xlim/30
  if (is.null(x_lab)) x_lab <- var
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])
  df$var <- df[,var]
  by_colors <- CreateColorsByAny(by=by, df=df)
  ## hist plot
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    xlab(x_lab) + ylab("Density") + labs(title = title) +
    scale_x_continuous(limits=c(0,xlim)) +
  geom_bar(aes(y = ..density.., fill=by), color="black", binwidth=bin_width)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  ## For "hold_axes"
  if (hold_axes == TRUE) {
  build <- ggplot_build(g)
  g <- g + coord_cartesian(ylim= c(min(build$panel$ranges[[1]]$y.range),
                                   max(build$panel$ranges[[1]]$y.range)))
  }
  ## Create and plot 'pars' (if fit_pareto = TRUE)
  if(fit_pareto == TRUE){
  if(is.null(pars)) pars <- FitParetoParsToData(df, var="var", by="by")
  if(is.null(by)){
    pars$by <- "all"
  } else {
    pars[,length(pars)+1] <- pars$by
    colnames(pars)[length(pars)] <- by
  }
  dens <- FitParetoDensityToArray(df=df, var=var, pars=pars, xlim=xlim)
  if(is.null(by)){
    dens$by <- "all"
  } else {
    dens[,length(dens)+1] <- dens$by
    colnames(dens)[length(dens)] <- by
  }
  g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
  }

  ## Set axes for the rest of the plots
  build <- ggplot_build(g)
  xmax <- max(build$panel$ranges[[1]]$x.range)
  ymax <- max(build$panel$ranges[[1]]$y.range)

  ## Create and plot 'quantiles'
  probs = c(0.5, 0.75, 0.95)
  quantiles <- ddply(df, "by", function(df) quantile(df$va,probs=probs,
    na.rm=TRUE))
  colnames(quantiles)[2:(length(probs)+1)] <- paste("q",probs, sep="")
  quantiles$xmax <- xmax  # for geom_text
  quantiles$ymax <- ymax  # for geom_text
  if (emp_labels == TRUE) {
  g <- g + geom_text(data = quantiles,
      aes(x = c(q0.5+(xmax*0.02), q0.75+(xmax*.02), q0.95+(xmax*.02)),
        y = c(ymax*.8, ymax*.6, ymax*.4),
        label = c(paste("Median:","\n",as.integer(q0.5)),
        paste("75th:","\n",as.integer(q0.75)),
        paste("95th:","\n",as.integer(q0.95)))),
      color = c("black", "gray20", "gray30"), hjust=0, vjust=1)
  }
  if (emp_lines == TRUE) {
  g <- g + geom_vline(data = quantiles, aes(xintercept = c(q0.5, q0.75, q0.95)),
      color = c("black", "gray20", "gray30"), size = 1,
      linetype = c("longdash", "dashed", "dotted"))
  }
  ## Plot 'pars' quantile lines and parameter value text
  if(fit_pareto == TRUE) {
  pars$xmax <- xmax  # for geom_text
  pars$ymax <- ymax  # for geom_text
  g <- g + geom_text(data = pars, aes(x = xmax*.95,y = ymax*.95,
    label = paste("Location: ", signif(location,3),"\n",
    "Scale: ", signif(scale,3),"\n",
    "Shape: ", signif(shape,3)), sep=""),
    color = fit_color, hjust=1, vjust=1)
  if (lines == TRUE) {
  g <- g + geom_vline(data = pars, aes(xintercept =
    c(VGAM::qgpd(.5, location=location, scale=scale, shape=shape),
      VGAM::qgpd(.75, location=location, scale=scale, shape=shape),
      VGAM::qgpd(.95, location=location, scale=scale, shape=shape))),
    linetype=c("longdash", "dashed", "dotted"), colour=fit_color, size=1)
  }
  if (labels == TRUE) {
  g <- g + geom_text(data = pars,
    aes(x = c((VGAM::qgpd(.5, location=location,scale=scale,shape=shape) +
      (xmax*0.02)), (VGAM::qgpd(.75, location=location, scale=scale,
      shape=shape) + (xmax*0.02)), (VGAM::qgpd(.95, location=location,
      scale=scale, shape=shape)) + (xmax*0.02)), y = c(ymax*.9, ymax*.7,
      ymax*.5),
    label = c(paste("Median:","\n",as.integer(VGAM::qgpd(.5, location=location,
      scale=scale, shape=shape))),
    paste("75th:","\n",as.integer(VGAM::qgpd(.75, location=location,
      scale=scale, shape=shape))),
    paste("95th:","\n",as.integer(VGAM::qgpd(.95, location=location,
      scale=scale, shape=shape))))),
    color= fit_color,hjust=0, vjust=1)
  }
  }
  g
}

#' Plots a histogram with a wrapped Cauchy distribution
#'
#' Plots a histogram of the data with an option to fit a wrapped Cauchy
#'   distribution
#'
#' Usage: PlotHistAndWrappedCauchy(df, var, by, pars, bin_width,
#'   fit_wrapped_cauchy, fit_color, x_lab, title)
#'
#' @param df Dataframe of data
#' @param var Column name of variable to fit wrapped Cauchy distribution
#' @param by Column name used to subset data, optional. Default is NULL.
#' @param pars A set of parameters used in place of pars generated within the
#'   functions, optional. Default is NULL.
#' @param bin_width Numeric, bin size. Default is 15 degrees or 2*pi/24 radians.
#' @param fit_wrapped_cauchy Logical, whether or not to fit and show wrapped
#'   Cauchy distribution. Default is TRUE.
#' @param fit_color String, color used for Pareto fit line, quantiles, and
#'   parameter value text. Default is "black".
#' @param x_lab String, name for x-axis. Default is 'var'.
#' @param title String, name for title.
#'
#' @return A plot of the data with a fitted wrapped Cauchy distribution.
#' @export
#'
#' @details Automatically adjusts plot for degrees or radians input, but all
#'   parameter estimates are based on radians.
#'
PlotHistAndWrappedCauchy<- function(df,
                                    var,
                                    by = NULL,
                                    pars = NULL,
                                    bin_width = NULL,
                                    fit_wrapped_cauchy = TRUE,
                                    fit_color = "black",
                                    x_lab = NULL,
                                    title = NULL) {
  ifelse(max(df[,var], na.rm=TRUE) <= 2*pi, radian <- TRUE, radian <- FALSE)
  if (is.null(x_lab)) x_lab <- var
  if (is.null(bin_width) & radian == TRUE) bin_width = (2*pi)/24
  if (radian == FALSE) {
    ifelse (is.null(bin_width), bin_width <- 15*(pi/180), bin_width <-
      bin_width*(pi/180))
  }
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])
  ifelse(radian == TRUE, df$var <- df[,var], df$var <- df[,var]*(pi/180))
  by_colors <- CreateColorsByAny(by=by, df=df)
  breaks <- seq(0, (2*pi), by=((2*pi)/12))
  breaks <- breaks[-length(breaks)]
  minor_breaks <- seq(0, 2*pi, by=bin_width)
  limits <- c(0, 2*pi)
  labels <- round(breaks, 2)
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    theme(axis.text.x = element_text(colour="grey20",size=12))+
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(vjust=3)) +
    xlab(x_lab) + ylab("Density") + labs(title = title)
  g <- g + geom_bar(aes(y = ..density.., fill=by), color="black",
    binwidth=bin_width) + coord_polar(start=(1.5*pi), direction = -1) +
    scale_x_continuous(limits=limits, breaks=breaks, minor_breaks=minor_breaks,
      labels = labels) + scale_y_continuous(labels = NULL)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  g
  if (fit_wrapped_cauchy == TRUE) {
    if (is.null(pars)) {
      pars <- FitWrappedCauchyParsToData(df=df, var="var", by="by")
      if(is.null(by)){
        pars$by <- "all"
      } else {
        pars[,length(pars)+1] <- pars$by
        colnames(pars)[length(pars)] <- by
      }
    }
    dens <- FitWrappedCauchyDensityToArray(var=var, pars=pars)
    if(is.null(by)){
      dens$by <- "all"
    } else {
      dens[,length(dens)+1] <- dens$by
      colnames(dens)[length(dens)] <- by
    }
    g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
    build <- ggplot_build(g)
    pars$xmax <- max(build$panel$ranges[[1]]$theta.range)  # for geom_text
    pars$ymax <- max(build$panel$ranges[[1]]$r.range)  # for geom_text
    g <- g + geom_text(data = pars, aes(y = ymax, x = xmax*.125,
      label = paste("mu:", signif(mu,3),"\n","rho:", signif(rho,3)), sep=""),
      size=4.5, color = fit_color, hjust=0, vjust=-1.5)
  }
  g
}

#' Plot histogram and wrapped Normal distribution
#'
#' Plots a histogram of the data with an option to fit a wrapped normal
#'   distribution.
#'
#' @usage PlotHistAndWrappedNormal (df, var, by, pars, bin_width,
#'   fit_wrapped_normal, fit_color, x_lab, title)
#'
#' @param df Dataframe of data
#' @param var Column name of variable to fit wrapped Cauchy distribution
#' @param by Column name used to subset data, optional. Default is NULL.
#' @param pars A set of parameters used in place of pars generated within the
#'   functions, optional. Default is NULL.
#' @param bin_width Numeric, bin size. Default is 15 degrees or 2*pi/24 radians.
#' @param fit_wrapped_normal Logical, whether or not to fit and show wrapped
#'   Normal distribution. Default is TRUE.
#' @param fit_color String, color used for Pareto fit line, quantiles, and
#'   parameter value text. Default is "black".
#' @param x_lab String, name for x-axis. Default is 'var'.
#' @param title String, name for title.
#'
#' @return A plot of the data with a fitted wrapped Normal distribution.
#' @export
#'
#' @details Automatically adjusts plot for degrees or radians input, but all
#'   parameter estimates are based on radians.
#'
PlotHistAndWrappedNormal <- function(df,
                                     var,
                                     by = NULL,
                                     pars = NULL,
                                     bin_width = NULL,
                                     fit_wrapped_normal = TRUE,
                                     fit_color = "black",
                                     x_lab = NULL,
                                     title = NULL) {
  ifelse(max(df[,var], na.rm=TRUE) <= 2*pi, radian <- TRUE, radian <- FALSE)
  if (is.null(x_lab)) x_lab <- var
  if (is.null(bin_width) & radian == TRUE) bin_width = (2*pi)/24
  if (radian == FALSE) {
    ifelse (is.null(bin_width), bin_width <- 15*(pi/180), bin_width <-
      bin_width*(pi/180))
  }
  ifelse(is.null(by), keep <- var, keep <- c(var,by))
  df <- subset(df, select = keep)
  ifelse(is.null(by),  df$by <- "all", df$by <- df[,by])
  ifelse(radian == TRUE, df$var <- df[,var], df$var <- df[,var]*(pi/180))
  by_colors <- suppressWarnings(CreateColorsByAny(by="by", df=df))
  breaks <- seq(0, (2*pi), by=((2*pi)/12))
  breaks <- breaks[-length(breaks)]
  minor_breaks <- seq(0, 2*pi, by=bin_width)
  limits <- c(0, 2*pi)
  if (radian == TRUE){
    labels <- round(breaks, 2)
  } else {
    labels <- round(breaks*(180/pi), 2)
  }
  g <- ggplot(df, aes(x=var)) +
    scale_fill_manual(values=by_colors) +
    theme(legend.position="none") +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    theme(axis.text.x = element_text(colour="grey20",size=12))+
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(vjust=3)) +
    xlab(x_lab) + ylab("Density") + labs(title = title, hjust=1)
  g <- g + geom_bar(aes(y = ..density.., fill=by), color="black",
    binwidth=bin_width) + coord_polar(start=(1.5*pi), direction = -1) +
    scale_x_continuous(limits=limits, breaks=breaks, minor_breaks=minor_breaks,
      labels = labels) + scale_y_continuous(labels = NULL)
  if (!is.null(by)) g <- g + facet_wrap( ~ by)
  if (fit_wrapped_normal == TRUE) {
    if (is.null(pars)) {
      pars <- FitWrappedNormalParsToData(df=df, var="var", by="by")
      if(is.null(by)){
        pars$by <- "all"
      } else {
        pars[,length(pars)+1] <- pars$by
        colnames(pars)[length(pars)] <- by
      }
    }
    dens <- suppressWarnings(FitWrappedNormalDensityToArray(var=var, pars=pars))
    if(is.null(by)){
      dens$by <- "all"
    } else {
      dens[,length(dens)+1] <- dens$by
      colnames(dens)[length(dens)] <- by
    }
    g <- g + geom_line(data=dens, aes(x=x, y=y), color=fit_color, size=1)
    build <- ggplot_build(g)
    pars$xmax <- max(build$panel$ranges[[1]]$theta.range)  # for geom_text
    pars$ymax <- max(build$panel$ranges[[1]]$r.range)  # for geom_text
    g <- g + geom_text(data = pars, aes(y = ymax, x = xmax*.125,
      label = paste("mu:", signif(mu,3),"\n","rho:", signif(rho,3)), sep=""),
      size=4.5, color = fit_color, hjust=0, vjust=-1.5)
  }
  g
}

#' Plot Logistic probability distribution function
#'
#' Plot Logistic Probability Distribution Function using logistic function with
#'   location and scale parameters
#'
#' @usage PlotLogisticPDF(location, scale, range, col)
#'
#' @param location Numeric, location parameter
#' @param scale Numeric, location parameter
#' @param range Vector for x-scale range, default is (0,50)
#' @param col String, color palette, e.g. "jet2.col(100)"
#'
#' @return Plot of probability distribution
#'
#' @import ggplot2
#' @export
#'
PlotLogisticPDF <- function(location,
                            scale,
                            range = c(0,50),
                            col = "darkgreen"){
  ifelse (length(range(range)) <= 100, length_x <- 100, length <- length(range))
  x <- seq(range[1], range[2], length=length_x)
  y <- exp(location + scale*x)/(1 + exp(location + scale*x))
  df <- as.data.frame(cbind(x, y))
  main <- paste0("Logistic (location = ", round(location, 5), ", scale = ",
    round(scale, 5), ")")
  g <- ggplot(df, aes(x, y)) +
    geom_line(aes(color=y), size=1.5) +
    scale_colour_gradientn(colours = col) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x='Value', y='Probability Density', title=main)
  g
}

#' Plot logistic probability density function
#'
#' Plot logistic probability density function using a logistic function with
#'   inflection and scale parameters
#'
#' @usage PlotLogisticInflectionPDF(inflection, scale, range, col)
#'
#' @param inflection Numeric, inflection parameter.
#' @param scale Numeric, scale parameter.
#' @param range Vector for x-scale range, default is (0,50).
#' @param col String, color palette, e.g. "jet2.col(100)".
#'
#' @import ggplot2
#'
#' @return Plot of probability distribution
#' @export
#'
PlotLogisticInflectionPDF <- function(inflection,
                                      scale,
                                      range = c(0,50),
                                      col = "darkgreen"){
  ifelse (length(range(range)) <= 100, length_x <- 100, length <- length(range))
  x <- seq(range[1], range[2], length=length_x)
  y <- LogisticByInflection(x, inflection, scale)
  df <- as.data.frame(cbind(x, y))
  main <- paste0("Logistic (inflection = ", round(inflection, 5), ", scale = ",
    round(scale, 5), ")")
  g <- ggplot(df, aes(x, y)) +
    geom_line(aes(color=y), size=1.5) +
    scale_colour_gradientn(colours = col) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x='Value', y='Probability Density', title=main)
  g
}

#' Plot Pareto probability density function
#'
#' Plot Pareto probability density function using the location, scale, and
#'   shape parameters
#'
#' @param location Numeric, location parameter, default is 1
#' @param scale Numeric, scale parameter, default is 1.
#' @param shape Numeric, shape parameter, default is 1.
#' @param xlim Numeric, sets plot x lim max, default is 99th percentile.
#' @param color
#'
#' @return Plot of Parteo probability density function
#' @export
#'
#' @details x-axis is set to have 500 values. Commented out section: 75th and
#'   95th percentiles are shown with dashed red lines.
PlotParetoPDF <- function(location = 1,
                          scale = 1,
                          shape = 1,
                          xlim = NULL,
                          color = "darkgreen") {
  if (is.null(xlim)){
    x <- seq(location, VGAM::qgpd(.99, location, scale, shape), length=501)
  } else {
    x <- seq(location, xlim, length=501)
  }
  x <- x[2:501]
  y <- VGAM::dgpd(x, location=location, scale=scale, shape=shape)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Pareto Distribution", "\n", "(location = ", signif(location,3),
    ", scale = ", signif(scale,3),
    ", shape = ", signif(shape,3), ")",  sep="")
  g <- ggplot(df,aes(x,y)) +
    geom_line(color=color, size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    scale_x_continuous (limits= c(0, max(x)))
  g + labs(x = 'Value', y= 'Probability Density', title=main)
#  geom_vline(xintercept = VGAM::qgpd(.75, location=location, scale=scale,
#    shape=shape), size=1, colour="red20", linetype = "longdash") +
#  geom_vline(xintercept = VGAM::qgpd(.95, location=location, scale=scale,
#    shape=shape), size=1, colour="red80", linetype = "longdash")
}

#' Plot Beta probability density functions
#'
#' Plot Beta probability density functions from list of parameter values
#'
#' @usage PlotParsBeta(pars, length)
#'
#' @param pars List, parameter values.
#' @param length Integer, number of x values used to plot lines, default is 100.
#'
#' @return Plot of BetaPDF with the parameter values of parameter list
#' @export
#'
#' @details Automatically adds 1 to the length, since 0 is included
PlotParsBeta <- function(pars,
                         length=100){
  behavior_colors <- CreateColorsByBehavior(output=TRUE)
  male_pars <- pars[["male"]]
  female_pars <- pars[["female"]]
  length = length + 1
  x <- seq(0, 1, length=length)
  behavior <- rep(NA, length)
  m <- rep("male",length)
  f <- rep("female", length)
  male <- data.frame(sex=m, x)
  female <- data.frame(sex=f, x)
  male$cruise <- dbeta(x, male_pars[["cruise_pars"]][["shape1"]],
                          male_pars[["cruise_pars"]][["shape2"]])
  male$forage <- dbeta(x, male_pars[["forage_pars"]][["shape1"]],
                          male_pars[["forage_pars"]][["shape2"]])
  male$nest <- dbeta(x, male_pars[["nest_pars"]][["shape1"]],
                        male_pars[["nest_pars"]][["shape2"]])
  male$loaf <- dbeta(x, male_pars[["loaf_pars"]][["shape1"]],
                        male_pars[["loaf_pars"]][["shape2"]])
  male$roost <- dbeta(x, male_pars[["roost_pars"]][["shape1"]],
                         male_pars[["roost_pars"]][["shape2"]])
  male$territorial <- dbeta(x,male_pars[["territorial_pars"]][["shape1"]],
                              male_pars[["territorial_pars"]][["shape2"]])
  female$cruise <- dbeta(x, female_pars[["cruise_pars"]][["shape1"]],
                            female_pars[["cruise_pars"]][["shape2"]])
  female$forage <- dbeta(x, female_pars[["forage_pars"]][["shape1"]],
                            female_pars[["forage_pars"]][["shape2"]])
  female$nest <- dbeta(x, female_pars[["nest_pars"]][["shape1"]],
                          female_pars[["nest_pars"]][["shape2"]])
  female$loaf <- dbeta(x, female_pars[["loaf_pars"]][["shape1"]],
                          female_pars[["loaf_pars"]][["shape2"]])
  female$roost <- dbeta(x, female_pars[["roost_pars"]][["shape1"]],
                           female_pars[["roost_pars"]][["shape2"]])
  female$territorial <- dbeta(x,female_pars[["territorial_pars"]][["shape1"]],
                                female_pars[["territorial_pars"]][["shape2"]])
  df <- as.data.frame(rbind(male, female))
  melted <- melt(df, c("sex","x"), 3:length(df), variable.name="behavior")
  ggplot(melted, aes(x = x, y=value, group= behavior, color=behavior)) +
    facet_grid(~ sex) + theme(panel.margin=unit(1, "lines")) +
    scale_y_continuous(limits=c(0, 5)) +
    scale_color_manual(values=behavior_colors) +
    geom_line(stat="identity", size=1.5) +
    scale_x_continuous(breaks=seq(0,1,.1)) +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) + labs(x="Daily Period",
    y="Probability Density", title="Daily Behavior Probability Distributions")
}

#' Plot two logistic probability density functions simultaneously
#'
#' Plot logistic probability distribution function using location and scale
#'   parameters.
#'
#' @usage PlotTwoLogisticPDF(location_1, scale_1, location_2, scale_2, x_label,
#'   x_max)
#'
#' @param location_1 Numeric, location parameter for distribution 1.
#' @param scale_1 Numeric, scale parameter for distribution 1.
#' @param location_2 Numeric, location parameter for distribution 2.
#' @param scale_2 Numeric, scale parameter for distribution 2.
#' @param x_label String, label for x-axis, default if "X Value".
#' @param x_max Numeric, maximum value on x scale, default is 1000.
#'
#' @return Plot of probability distribution
#' @export
#'
PlotTwoLogisticCDF <- function(location_1,
                               scale_1,
                               location_2,
                               scale_2,
                               x_label = "X Value",
                               x_max = 1000) {
  x <- seq(0, x_max, length=1000)
  y <- exp(location_1 + scale_1*x)/(1 + exp(location_1 + scale_1*x))
  y2 <- exp(location_2 + scale_2*x)/(1 + exp(location_2 + scale_2*x))
  df <- as.data.frame(cbind(x, y, y2))
  df$x_max <- x_max
  g <- ggplot(df) +
    geom_line(aes(x, y), colour="blue3", size=1.5) +
    geom_line(aes(x, y2), colour="red3", size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x=x_label, y='Probability Density') +
    annotate("text", x=x_max*.05, y=.5, label=paste0("Location: ",
      round(location_1, 5),"\n", "Scale: ", round(scale_1, 5)), size=6,
      color="blue3", hjust=0, vjust=1) +
    annotate("text", x=x_max*.95, y=.5, label=paste0("Location: ",
      round(location_2, 2),"\n", "Scale: ", round(scale_2, 5)), size=6,
      color="red3", hjust=1, vjust=1)
  g
}

#' Plot two logistic probability density function
#'
#' Plot two logistic probability density function (using inflection and scale
#'   parameters)
#'
#' @param inflection_1 Numeric, inflection parameter for distribution 1.
#' @param scale_1 Numeric, scale parameter for distribution 1.
#' @param inflection_2 Numeric, inflection parameter for distribution 2.
#' @param scale_2 Numeric, scale parameter for distribution 2.
#' @param x_label String, label for x-axis, default if "X Value".
#' @param log1_range Vector, log1 value range. Default = c(0, 45).
#' @param log2_range Vector, log2 value range. Default = c(0, 45).
#' @param xlab_range Vector, x value range. Default = c(0, 45).
#'
#' @return Plot of probability distributions.
#' @export
#'
PlotTwoLogisticInflectionPDF <- function(inflection_1,
                                         scale_1,
                                         inflection_2,
                                         scale_2,
                                         x_label = "X Value",
                                         log1_range = c(0, 45),
                                         log2_range = c(0, 45),
                                         xlab_range = c(0, 45)) {
  log1_seq <- seq(log1_range[1], log1_range[2], length=1000)
  log2_seq <- seq(log2_range[1], log2_range[2], length=1000)
  x <- seq(xlab_range[1], xlab_range[2], length=1000)
  x_max <- xlab_range[2]
  y <- LogisticByInflection(log1_seq, inflection_1, scale_1)
  y2 <- LogisticByInflection(log2_seq, inflection_2, scale_2)
  log1 <- cbind.data.frame(log1_seq, y)
  log2 <- cbind.data.frame(log2_seq, y2)
  df <- as.data.frame(cbind(x, y, y2))
  df$x_max <- x_max
  g <- ggplot() +
    geom_line(data = log1, aes(log1_seq, y), colour="slateblue4", size=1.5) +
    geom_line(data = log2, aes(log2_seq, y2), colour="saddlebrown", size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    labs(x=x_label, y='Probability Density') # +
#    annotate("text", x=x_max*.05, y=.5, label=paste0("Inflection: ",
#      signif(inflection_1, 3),"\n", "Scale: ", signif(scale_1, 4)), size=6,
#      color="blue3", hjust=0, vjust=1) +
#    annotate("text", x=x_max*.95, y=.5, label=paste0("Inflection: ",
#      signif(inflection_2, 3),"\n", "Scale: ", signif(scale_2, 4)), size=6,
#      color="red3", hjust=1, vjust=1)
  g
}

#' Plot Weibull cumulative distribution function
#'
#' Plot Weibull cumulative distribution function using shape and scale
#'   parameters
#'
#' @usage PlotWeibullCDF(shape, scale, max_x)
#'
#' @param shape Numeric, Weibull shape, default is 1.
#' @param scale Numeric, Weibull scale, default is 1.
#' @param max_x Numeric, maximum value on x scale, default is 250.
#'
#' @return Plot of probability distribution
#' @export
#'
#' @import ggplot2
#'
PlotWeibullCDF <- function(shape = 1,
                           scale = 1,
                           max_x = 250){
  x <- seq(0, max_x, length=250)
  y <- pweibull(x, shape=shape, scale=scale)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Weibull Distribution (shape = ", shape, ", scale = ",
    scale, ")",  sep="")
  g <- ggplot(df, aes(x, y)) +
    geom_line(colour="dark green", size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x='Value', y='Cumulative Probability Density', title=main)
}

#' Plot Weibull cumulative distribution function
#'
#' Plot Weibull cumulative distribution function using shape and scale
#'   parameters
#'
#' @usage PlotWeibullCDF(shape, scale, max_x)
#'
#' @param shape Numeric, Weibull shape, default is 1.
#' @param scale Numeric, Weibull scale, default is 1.
#' @param max_x Numeric, maximum value on x scale, default is 250.
#'
#' @return Plot of probability distribution
#' @export
#' @import ggplot2
#'
PlotWeibullPDF <- function(shape = 1,
                           scale = 1,
                           max_x = 250){
  x <- seq(0, max_x, length=251)
  y <- dweibull(x, shape=shape, scale=scale)
  df <- as.data.frame(cbind(x, y))
  main <- paste("Weibull Distribution (shape = ", shape, ", scale = ",
    scale, ")", sep="")
  g <- ggplot(df, aes(x, y)) +
    geom_line(colour="dark green", size=1.5) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black"))
  g + labs(x='Value', y='Probability Density', title=main)
}

#' Summarizes a 'fitdist' object
#'
#' Summarizes a 'fitdist' object (from 'fitdistrplus' package), made to be used
#'   as a table (e.g. xtable)
#'
#' @usage SummarizeFistdist(fits)
#'
#' @param fits 'fitdist' object.
#'
#' @return Dataframe
#' @export
#'
SummarizeFitdist <- function(fits) {
  fits <- fits
  df <- as.data.frame(cbind(c("parameter", "estimate", "sd", "loglik",
    "aic", "bic", "method", "ddistname", "pdistname", "qdistname")))
  colnames(df)[1] <- NA
  for (i in 1:length(fits)){
    start_col <- length(df) + 1
    n_pars <- length(fits[[i]]$estimate)
    for (j in 1:n_pars) df <- cbind(df, NA)
    colnames(df)[start_col] <- names(fits[i])
    for (k in 1:n_pars) {
      df[1, start_col + k - 1] <- names(fits[[i]]$estimate)[k]
      df[2, start_col + k - 1] <- signif(fits[[i]]$estimate[k], 4)
      df[3, start_col + k - 1] <- signif(fits[[i]]$sd[k], 4)
    }
    df[4, start_col] <- format(fits[[i]]$loglik, nsmall=2)
    df[5, start_col] <- format(fits[[i]]$aic, nsmall=2)
    df[6, start_col] <- format(fits[[i]]$bic, nsmall=2)
    df[7, start_col] <- fits[[i]]$method
    df[8, start_col] <- paste0("d", fits[[i]]$distname)
    df[9, start_col] <- paste0("p", fits[[i]]$distname)
    df[10, start_col] <- paste0("q", fits[[i]]$distname)
  }
  colnames(df)[1] <- ""
  colnames(df)[which(names(df) == "NA")] <- ""
  return(df)
}


#' Summarize data and standard error
#'
#' Summarizes data with count, mean, standard deviation, standard error of the
#'    mean, and confidence interval
#'
#' @usage SummarizeSE(df, var, groups, na_rm, conf_int)
#'
#' @param data Dataframe
#' @param var Column name that contains the variable to be summarized
#' @param by Column name of vector containing grouping variables, Default is
#'   NULL.
#' @param na_rm Logical that indicates whether to ignore NAs.
#' @param conf_interval Numeric, the percent range of the confidence interval.
#'   Default is .95)
#'
#' @return Dataframe
#' @export
#'
#' @details Based on example function from "Cookbook for R" website.
SummarizeSE <- function(data,
                        var,
                        by = NULL,
                        na_rm = FALSE,
                        conf_interval = .95) {
  LengthNaRm <- function (x, na.rm=FALSE) {
    if (na_rm == TRUE){
      sum(!is.na(x))
    } else {
      length(x)
    }
  }
  ifelse(is.null(by), data$by <- "all", data$by <- data[,by])
  data2 <- plyr::ddply(data, by, .fun = function(x, var) {c(n =
    LengthNaRm(x[[var]], na.rm=na_rm), mean = mean(x[[var]], na.rm=na_rm),
    variance = var(x[[var]], na.rm=na_rm), sd = sd(x[[var]], na.rm=na_rm))},
    var)
  data2 <- rename(data2, c("mean" = var))
  data2$se <- data2$sd / sqrt(data2$n)  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ci <- qt(conf_interval/2 + .5, data2$n-1)
  data2$ci <- data2$se * ci
  return(data2)
}
