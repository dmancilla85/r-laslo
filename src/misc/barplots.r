########################################
########################################
#
# Building barplots in R with
# error bars
#
########################################
########################################

########################################
# Constructing our data set
########################################

# Aggregate our data from the mtcars dataset
# grouping by gears and cylinders
#
# Our custom function will return the mean
# mileage and its standard deviation for each
# combination of gears and cylinders, as well
# as the number of observations for that
# combination of gears and cylinders
myData <- aggregate(mtcars$mpg,
  by = list(cyl = mtcars$cyl, gears = mtcars$gear),
  FUN = function(x) {
    c(
      mean = mean(x), sd = sd(x),
      n = length(x)
    )
  }
)

# Turn the result of the last function
# into something useable
myData <- do.call(data.frame, myData)

# Calculate the standard error
myData$se <- myData$x.sd / sqrt(myData$x.n)

# Rename the columns for ease of use
colnames(myData) <- c("cyl", "gears", "mean", "sd", "n", "se")

# Construct a group column for plot legends
myData$names <- c(paste(
  myData$cyl, "cyl /",
  myData$gears, " gear"
))


########################################
# Construct a basic non-grouped
# barplot for each set of observations
########################################

# Specify plot margins
par(mar = c(5, 6, 4, 5) + 0.1)

# Set custom upper limits for the plot
# so that our error bars don't run into
# the title
plotTop <- max(myData$mean) +
  myData[myData$mean == max(myData$mean), 6] * 3

# Construct a barplot with heights from
# our data
# Specify names of each column with names.arg
# beside = true for non-stacked bars
# las = 2 for all horizontal axis labels
barCenters <- barplot(
  height = myData$mean,
  names.arg = myData$names,
  beside = true, las = 2,
  ylim = c(0, plotTop),
  cex.names = 0.75, xaxt = "n",
  main = "Mileage by No. Cylinders and No. Gears",
  ylab = "Miles per Gallon",
  border = "black", axes = TRUE
)

# Specify the groupings. We use srt = 45 for a
# 45 degree string rotation
text(
  x = barCenters, y = par("usr")[3] - 1, srt = 45,
  adj = 1, labels = myData$names, xpd = TRUE
)

# Add in line segments and caps for error bars
# set to plus/minus 2 standard errors from the mean
segments(barCenters, myData$mean - myData$se * 2, barCenters,
  myData$mean + myData$se * 2,
  lwd = 1.5
)

arrows(barCenters, myData$mean - myData$se * 2, barCenters,
  myData$mean + myData$se * 2,
  lwd = 1.5, angle = 90,
  code = 3, length = 0.05
)


########################################
# Group our bar graphs and prettify
# with some custon styling
########################################

# Turn our means and standard errors from
# vectors to matrices to group them
tabbedMeans <- tapply(
  myData$mean, list(myData$cyl, myData$gears),
  function(x) c(x = x)
)

tabbedSE <- tapply(
  myData$se, list(myData$cyl, myData$gears),
  function(x) c(x = x)
)

# Make sure our y axis goes high enough
plotTop <- max(myData$mean) +
  myData[myData$mean == max(myData$mean), 6] * 3

# Construct our barplot using our means matrix
# Specify some custom colors
# Play around with legend positioning
barCenters <- barplot(
  height = tabbedMeans,
  beside = TRUE, las = 1,
  ylim = c(0, plotTop),
  cex.names = 1.15,
  cex.axis = 1,
  cex.lab = 1.15,
  col = c("steelblue", "steelblue2", "steelblue4"),
  main = "Mileage by No. Cylinders and No. Gears",
  ylab = "Miles per Gallon",
  xlab = "Number of  Gears",
  border = "black", axes = TRUE,
  legend.text = TRUE,
  args.legend = list(
    title = "No. Cylinders",
    x = 15,
    y = 30,
    bty = "y",
    cex = .7
  )
)

# Add in line segments and caps for error bars
segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
  tabbedMeans + tabbedSE * 2,
  lwd = 1.5
)

arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
  tabbedMeans + tabbedSE * 2,
  lwd = 1.5, angle = 90,
  code = 3, length = 0.05
)


########################################
# Recreate ugly non-grouped barplot
# using ggplot2
########################################

# load in ggplot2 library
library(ggplot2)

# Set spacing between bars
dodge <- position_dodge(width = 0.9)

# Set y limits
limits <- aes(
  ymax = myData$mean + myData$se,
  ymin = myData$mean - myData$se
)

# Construct ggplot object
p <- ggplot(data = myData, aes(x = names, y = mean, fill = names))

# Add in bars, error bars, and fiddle with
# axis names and ticks and titles
p + geom_bar(stat = "identity", position = dodge) +
  geom_errorbar(limits, position = dodge, width = 0.25) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(x = "No. Cylinders and Gears", y = "Miles Per Gallon") +
  ggtitle("Mileage by No. Cylinders\nand No. Gears")


########################################
# Group our bar graphs in ggplot2
########################################

# Construct ggplot object specifying cylinders
# and gears as factors to make sure they group
p <- ggplot(data = myData, aes(
  x = factor(cyl), y = mean,
  fill = factor(gears)
))

# Same stuff as last time: add in bars and
# error bars and titles and labels
p + geom_bar(
  stat = "identity",
  position = position_dodge(0.9)
) +
  geom_errorbar(limits,
    position = position_dodge(0.9),
    width = 0.25
  ) +
  labs(x = "No. Cylinders", y = "Miles Per Gallon") +
  ggtitle("Mileage by No. Cylinders\nand No. Gears") +
  scale_fill_discrete(name = "No. Gears")
