################################################################################
# This script creates csv of lake ids viable to run with lakeflow
# and sets up the saved filepaths to be compatible with job array
################################################################################
# Load module
library(data.table)

# Load the csvs of viable locations
files <- Sys.glob("in/viable/viable_locations*.csv") #note: may need to change this based on file names!
lakes <- rbindlist(lapply(files, fread))
colnames(lakes) <- "lake"
print(nrow(lakes))

# Break into 20 separate chunks
n <- 30
breaks <- round(seq(from=1, to=nrow(lakes), length.out=(n+1)))
for(i in 1:n){
    dt <- lakes[breaks[i]:breaks[i+1],]
    fwrite(dt[,"lake"], paste0("in/viable/lakeflow", i, ".csv"))
}
