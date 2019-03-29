#!/usr/bin/env Rscript
USAGE <- "Usage: $ fitdist.R families.counts.tsv [outdir]"
IMG_WIDTH = 1280
IMG_HEIGHT = 960
FONT_SIZE = 18
# Get arguments.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1 || args[1] == '-h') {
  print(paste(USAGE), quote=FALSE)
  quit()
}
infile <- args[1]
if (length(args) >= 2) {
  outdir <- args[2]
} else {
  outdir <- '.'
}
library(fitdistrplus)
library(logspline)

save_plot <- function(data, plot_fxn, dir, filename) {
  png(paste(dir, filename, sep='/'), width=IMG_WIDTH, height=IMG_HEIGHT, pointsize=FONT_SIZE)
  if (plot_fxn == 'plot') {
    plot(data)
  } else if (plot_fxn == 'descdist') {
    descdist(data, discrete=FALSE)
  }
  dev.off()
}

data <- read.table(infile)[,1]
fakedata <- data[data!=1]
save_plot(data, 'descdist', outdir, 'cullen-frey.raw.png')
save_plot(fakedata, 'descdist', outdir, 'cullen-frey.fake.png')

datas <- list(data, fakedata)
data_names <- c('raw', 'fake')
titles <- c('raw', 'singleton-free')
for (dist_name in c('exp', 'gamma', 'weibull')) {
  for (i in c(1,2)) {
    fit <- fitdist(datas[[i]], dist_name)
    save_plot(fit, 'plot', outdir, paste('fit', data_names[i], dist_name, 'png', sep='.'))
    #mtext(dist_name, side=3, line=20, adj=-0.3)
  }
}
