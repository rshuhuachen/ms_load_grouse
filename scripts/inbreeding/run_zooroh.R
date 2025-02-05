### load packages
install.packages("RZooRoH", repos='http://cran.us.r-project.org')
library(stringr); library(RZooRoH)

args <- commandArgs(trailingOnly = TRUE)
file_gen <- args[[1]]
#samples <- str_remove(args[[2]], "^\\../..\\/")
file_out_roh <- args[[2]]

## load in data
data <- zoodata(genofile = file_gen, zformat = "gp", chrcol = 1, poscol = 3, supcol = 5)

### define model
model <- zoomodel(krates=c(2,20,200,2000,2000), err = 0.005, K = 5)

### run model
out <- zoorun(zoomodel = model, 
              zooin = data, 
              nT = 8, localhbd = TRUE)

print("model done running")

save(out, file = "/output/inbreeding/zooroh_out.RData")

library(dplyr)

### divide out rzooroh ###
hbd <- realized(out)
hbd$file_id <- out@sampleids

save(hbd, file = "output/inbreeding/hbd_froh.RData")


