library(ggplot2)
DC_report <- read.table("output/genotyping/DC_report.ldepth.mean", sep = "\t", header = T)

ggplot(DC_report, aes (x = MEAN_DEPTH)) + geom_histogram() + theme_classic()

summary(DC_report$MEAN_DEPTH)  

#filter minimum coverage = 0.3 * mean and max coverage = 2*mean
mean(DC_report$MEAN_DEPTH)/3 # = 10.5 = 10
2*mean(DC_report$MEAN_DEPTH) # = 63
