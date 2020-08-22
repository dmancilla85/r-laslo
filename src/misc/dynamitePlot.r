dynamitePlot <- function(height, error, names = NA, 
        significance = NA, ylim = c(0,maxLim), …){
    maxLim <- 1.1* max(mapply(sum, height, error))
    bp <- barplot(height, names.arg = names, ylim = ylim, ...)
    arrows(x0 = bp, y0 = height, y1 = height + error, angle = 90)
    text(x = bp, y = 0.2 + height + error, labels = significance)
}
Values <- c(1,2,5,4)
Errors <- c(0.25, 0.5, 0.33, 0.12)
Names <- paste("Trial", 1:4)
Sig <- c("a", "a", "b", "b") dynamitePlot(Values, Errors, names = Names, significance = Sig)