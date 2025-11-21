################################################################################
# This script creates a csv of the lake ids to download data for lakeflow
# and sets up the saved filepaths to be compatible with job array
################################################################################
# Load module
library(data.table)

# Load PLD dataset with all lakes where LakeFlow could theoretically be run 
updated_pld <- fread("in/SWORDv16_PLDv103_wo_ghost_rch.csv")
updated_pld$lake_id <-  as.character(updated_pld$lake_id)
updated_pld$continent <- substr(updated_pld$lake_id, 1,1)

# Create a folder to store the lake id csvs
dir.create("in/lakeids", showWarnings = FALSE)

# Break into 50 separate chunks
n <- 50
breaks <- round(seq(from=1, to=nrow(updated_pld), length.out=(n+1)))
for(i in 1:n){
    dt <- updated_pld[breaks[i]:breaks[i+1],]
    fwrite(dt[,"lake_id"], paste0("in/lakeids/lakeid", i, ".csv"))
}
