#' Plots selected color palettes
#'
#' Plots selected color palettes (e.g. R_pal, SAGA_pal) with indexes and names
#'
#' @usage PlotColorPalette(pal, sel)
#'
#' @param pal Color palette, not in quotes (e.g. R_pal, SAGA_pal)
#' @param sel Vector, layer(s) in color palette to be displayed
#'
#' @return Plot of color palette(s)
#' @export
#'
#' @details Plot of color palette(s)
#'
PlotColorPalette <- function (pal,
                              sel = 1:length(pal)){
  for (i in 1:length(pal)){
    names(pal)[i] <- paste("[", i, "] ", names(pal)[i], sep = "")  # index num
  }
  if (length(sel) > 10){
    sel <- sel[1:10]
  }
  par(mfrow = c(length(sel), 1), mar = c(1.5, 0.8, 1.5, 0.5))
  for (j in sel){
    plot(y = rep(1, length(pal[[j]])), x = 1:length(pal[[j]]),
    axes = FALSE, xlab = "", ylab = "", pch = 15, cex = 5, col = pal[[j]])
    mtext(names(pal)[j], cex = 1, line = .1, side = 3)
  }
}

#' Plots colors in pie chart
#'
#' Plots selected colors from palettes in a pie chart
#'
#' @param pal the color palette (e.g. R_pal, SAGA_pal)
#' @param names Logical, whether or not to show name assigned to each color.
#'   Default is TRUE.
#' @param radius Numeric, the pie is drawn centered in a square box whose sides
#'   range from -1 to 1. If the character strings labeling the slices are long
#'   it may be necessary to use a smaller radius. Default is 1.
#'
#' @return Plot of color pie chart
#' @export
#'
#' @details Must attach appropriate color palettes prior to running function.
#'
PlotColorPie <- function(pal,
                         names = TRUE,
                         radius = 1){
  plot.new()
  par(mfrow=c(1,1))
  if (names == TRUE){
    pie(rep(1,length(pal)), labels = names(pal), col = pal, radius = radius)
  } else {
    pie(rep(1,length(pal)), labels = NA, col = pal, radius = radius)
  }
}

#' Plot histogram by grouping variable
#'
#' Plots a histogram of a variable by id of a grouping variable
#'
#' @param df Dataframe of data
#' @param var Column name with data to create histogram
#' @param id Column name of unique identified. Default is "id".
#'
#' @return A histogram plot
#' @export
#'
PlotHistogramByVar <- function(df,
                               var,
                               id = "id") {
  df <- df
  sum_df <- SummarizeSE(df, var, id)
  id_colors <- CreateIDColors(output=TRUE)
  grid <- seq(min(df[,var], na.rm = TRUE), max(df[,var], na.rm = TRUE),
      length = 100)
  normaldens <- ddply(df, id, function(df) {data.frame(var = grid,
    density = dnorm(grid, mean(df[,var]), sd(df[,var])))})
  g <- ggplot(df, aes(x = var, fill=id)) + facet_wrap( ~ id)  +
    scale_fill_manual(values=id_colors) +
    theme(legend.position="none") +
    theme(plot.title=element_text(size=22)) +
    theme(text=element_text(size=20, colour="black")) +
    theme(axis.text=element_text(colour="black")) +
    xlab("Speed") + ylab("Density")
  g + geom_bar(aes(y = ..density.., fill=id), colour="black", binwidth = 2) +
#    geom_text(aes(label=paste("mean: ", signif(speed,3),
#      "\n","sd: ",signif(sd,3), sep="")), data=sum_df) +
    geom_line(aes(y = density), colour="black", size=1, data = normaldens)
}

#' Removes all objects in global environment except selected ones.
#'
#' Removes all objects in global environment except ones included in the
#'   argument in quotes, e.g., c("df1", "vec1").
#'
#' @param object Object(s) to keep in environment, e.g., c("df1", "vec1").
#'
#' @return New global environment
#' @export
#'
RemoveExcept <- function(object = object){
  if (length(setdiff(ls(pos = .GlobalEnv), object)) > 0) {
    rm(list=setdiff(ls(pos = .GlobalEnv), object), pos = .GlobalEnv)
  }
}

#' Replace text in files
#'
#' Replaces a text string within files
#'
#' @usage ReplaceFilesText(files, text, replace)
#'
#' @param files List of filenames
#' @param text String of text to search for in files
#' @param replace String for replacement in files
#'
#' @return Files with replaced text
#' @export
#'
ReplaceFilesText <- function(files,
                             text,
                             replace) {
  for(i in files){
    x <- readLines(i)
    y <- gsub(text, replace, x)
    cat(y, file=i, sep="\n")
  }
}

#' Save last printed ggplot
#'
#' A wrapper function for ggsave()
#'
#' @usage SaveGGPlot(filename, path, width, height, units, dpi, bg)
#'
#' @param filename String, file name/filename of plot
#' @param path String, path to save plot to (if you just want to set path and
#'   not filename). If NULL(default), reverts to working directory
#' @param width Numeric, width. Default is 10
#' @param height Numeric, height. Default is 7.5
#' @param units String, units for width and height when either one is
#'   explicitly specified (in, cm, or mm). Default is "in".
#' @param dpi Integer, dpi to use for raster graphics. Default is 300.
#' @param bg String, background color. Default is "white"
#'
#' @return Saves a jpeg file of the last displayed plot
#' @export
#'
#' @details Default output format is set for PowerPoint presentations
#'
SaveGGPlot <- function (filename,
                        path = NULL,
                        width = 10,
                        height = 7.5,
                        units = "in",
                        dpi=300,
                        bg = "white"){
  if(is.null(path)) path <- getwd()
  if (file_ext(filename) == "") filename <- paste0(filename,".png")
  ggplot2::ggsave(filename = filename, path = path, width=width, height=height,
    units=units, dpi=dpi, bg = bg)
}

#' Saves last plot
#'
#' A function for saving the last displayed plot as a png, jpeg, or pdf.
#'
#' @usage SavePlot(filename, path, width, height, units, dpi)
#'
#' @param filename String, file name/filename of plot.
#' @param path String, path to save plot to (if you just want to set path and
#'   not filename). If NULL(default), uses working directory.
#' @param width Numeric, width. Default is 10.
#' @param height Numeric, height. Default is 7.5.
#' @param units Numeric, units for image width and height when either one is
#'   explicitly specified (in, cm, or mm), default is "in".
#' @param dpi Integer, dpi to use for plotting. Default is 300.
#'
#' @return Saves a file of the last displayed plot.
#' @export
#'
#' @examples Default output format width and height are set for PowerPoint
#'   presentations
#'
SavePlot <- function (filename,
                      path = NULL,
                      width = 10,
                      height = 7.5,
                      units = "in",
                      dpi = 300){
  if (is.null(path)) path <- getwd()
  if (file_ext(filename) == "") filename <- paste0(filename,".png")
  filepath <- file.path(path, filename)
  if (file_ext(filename) == "png") {
    dev.copy(png, filename=filepath, width=width, height=height, units=units,
      pointsize=12, bg="white", res=dpi)
    dev.off()
  }
  if (file_ext(filename) == "jpeg" || file_ext(filename) == "jpg") {
    dev.copy(jpeg, filename=filepath, width=width, height=height, res=dpi,
      units=units, quality=100)
    dev.off()
  }
  if (file_ext(filename) == "pdf") {
    dev.copy2pdf(file=filepath, width=width, height=height, paper="special")
  }
}

#' Black theme for ggplot
#'
#' A function that creates a black background/white text theme for ggplot
#'
#' @param base_size Numeric, font base size. Default is 12.
#' @param base_family String, font name. Default is "".
#'
#' @return A ggplot theme
#'
#' @import ggplot2
#' @export
#'
#' @details Original script is from: \url{https://gist.github.com/jslefche}
#'
ThemeBlack = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),
      axis.text.x = element_text(size = base_size*0.8, color = "white",
        lineheight = 0.9),
      axis.text.y = element_text(size = base_size*0.8, color = "white",
        lineheight = 0.9),
      axis.ticks = element_line(color = "white", size  =  0.2),
      axis.title.x = element_text(size = base_size, color = "white",
        margin = margin(0, 10, 0, 0)),
      axis.title.y = element_text(size = base_size, color = "white",
        angle = 90, margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold",
        hjust = 0, color = "white"),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_rect(fill = NA, color = "white"),
      panel.grid.major = element_line(color = "grey35"),
      panel.grid.minor = element_line(color = "grey20"),
      panel.margin = unit(0.5, "lines"),
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",
        angle = -90),
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")
    )
}
