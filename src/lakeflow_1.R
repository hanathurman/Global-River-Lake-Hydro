################################################################################
# This script downloads and saves data to run lakeflow
################################################################################

################################################################################
# Load libraries
################################################################################
library(ncdf4)
library(foreign)
library(lubridate)
library(data.table)
library(dplyr)
library(sf)
library(raster)
library(rstudioapi)
library(jsonlite)
library(httr)
library(BBmisc)
library(reticulate)
library(geoBAMr)
library(future)
library(future.apply)
library(optparse)
library(paws)
library(logger)
'%!in%' <- function(x,y)!('%in%'(x,y))

# Add in which python version to use with reticulate
use_python("/usr/local/bin/python3.12")

# Example Deployment using docker
# docker run --mount type=bind,source=C:/Users/kmcquil/Documents/LakeFlow_Confluence_Dev,target=/app lakeflow Rscript src/lakeflow_1.R "in/lakeids/lakeid1.csv" 1
# docker run -v /mnt/lakeflow:/data/input lakeflow_input -c /data/input/test/lakeid1.csv -w 1 -i /data/input

#use_virtualenv("r-reticulate")

################################################################################
# Constants

SWORD_VERSION = "16"

################################################################################
# Set args
option_list <- list(
    make_option(c("-c", "--input_file"), type = "character", default = NULL, help = "filepath to csv with lake ids to download data, point it to reaches_of_interest.json to process lakes associated with those reaches"),
    make_option(c("-w", "--workers"), type = "integer", default = NULL, help = "number of workers to use to download swot data"),
    make_option(c("-i", "--indir"), type = "character", default = NULL , help = "directory with input files"),
    make_option(c("-p", "--prefix"), type = "character", default = "", help = "prefix for hydrocron api-key storage"),
    ## I had an index argument, but the script is actually faster in seriel instead of calling hydrocron in parallel
    make_option(c("--index"), type = "integer", default = NULL , help = "Chooses what lake to process from input file, if -256 it uses array number")
  )
################################################################################

# Grab arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


# Number of workers used to download SWOT data
workers_ <- opts$workers

# Path that holds input data, was previously 'in'
indir <- opts$indir

# Prefix for hydrocron API in the parameter store
prefix <- opts$prefix

# Set index, use array if index is -256
index <- opts$index
if (index == -256){
  index <- strtoi(Sys.getenv("SLURM_ARRAY_TASK_ID"))
}

cat("SLURM_ARRAY_TASK_ID inside R:", Sys.getenv("SLURM_ARRAY_TASK_ID"), "\n")

################################################################################
# Load datasets
################################################################################

# Load PLD dataset
updated_pld = fread(file.path(indir,"/SWORDv16_PLDv103_wo_ghost_rch.csv"))
updated_pld$lake_id =  as.character(updated_pld$lake_id)
updated_pld$continent = substr(updated_pld$lake_id, 1,1)

# Load ET dataset
et = fread(file.path(indir, '/ancillary/et.csv'))
et$lake_id <- as.character(et$lake_id)

# Load tributary dataset
tributary = fread(file.path(indir,'/ancillary/tributaries.csv'))
tributary$lake_id <- as.character(tributary$lake_id)

# Load geoglows dataset
sword_geoglows = fread(file.path(indir,'/ancillary/sword_geoglows.csv'))
sword_geoglows$reach_id = as.character(sword_geoglows$reach_id)

# Create a folder to store the downloaded/processed datasets
dir.create(file.path(indir, "/clean"), showWarnings = FALSE)

################################################################################
# Define functions
################################################################################

get_connected_lake_ids <- function(reach_ids, pld) {
  # Ensure reach_ids are character (to match string fields)
  reach_ids <- as.character(reach_ids)
  
  # Function to check if any of the reach_ids are in a comma-separated field
  has_overlap <- function(field, reach_ids) {
    sapply(field, function(x) {
      ids <- unlist(strsplit(x, ","))
      any(ids %in% reach_ids)
    })
  }

  # Find rows where either upstream or downstream reach_id lists contain any of the reach_ids
  is_connected_up <- has_overlap(pld$U_reach_id, reach_ids)
  is_connected_dn <- has_overlap(pld$D_reach_id, reach_ids)

  # Combine the logical vectors and get lake_ids
  connected_lakes <- unique(pld$lake_id[is_connected_up | is_connected_dn])

  return(connected_lakes)
}

get_api_key <- function(prefix) {
    ssm <- paws::ssm()
    
    tryCatch({
      param_name <- paste0(prefix, "-hydrocron-key")
      response <- ssm$get_parameter(
        Name = param_name,
        WithDecryption = TRUE
      )
      api_key <- response$Parameter$Value
      log_info("Querying with Hydrocron API key.")
      return(api_key)
    }, error = function(e) {
      log_error(e$message)
      log_info("Not querying with Hydrocron API key.")
      return("")
    })
  }

pull_lake_data <- function(feature_id, api_key){
  website = paste0('https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?feature=PriorLake&feature_id=',feature_id, '&start_time=2023-01-01T00:00:00Z&end_time=2025-12-31T00:00:00Z&output=csv&fields=lake_id,time_str,wse,area_total,xovr_cal_q,partial_f,dark_frac,ice_clim_f,xtrk_dist,quality_f,partial_f')
  if (nzchar(api_key)) {
    # do something
    response = GET(website, add_headers("x-hydrocon-key" = api_key))
  } else {
    # do something else
    response = GET(website)
  } 
  # print(content(response))
  pull = content(response, as='parsed')$results
  if (!is.null(pull)) {
    data = try(read.csv(textConnection(pull$csv), sep=','))
  } else {
    return(NA)
  }

  if(is.error(data)){return(NA)}

  data$reach_id = feature_id
  return(data)
}

batch_download_SWOT_lakes <- function(obs_ids, api_key){
  print("batch downloading swot lakes")
  plan(multisession, workers = workers_)
  SWOT_data = future_lapply(unique(obs_ids),pull_lake_data, api_key = api_key)
  plan(sequential)
  print("finished batch downloading swot lakes")
  return(SWOT_data)
} 

pull_data <- function(feature_id, api_key){
  # print("pulling data")
  # Function to pull swot reach data using hydrocron
  website = paste0('https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries?feature=Reach&feature_id=',feature_id, '&start_time=2023-01-01T00:00:00Z&end_time=2025-12-31T00:00:00Z&output=csv&fields=reach_id,time_str,wse,width,slope,slope2,d_x_area,area_total,reach_q,p_width,xovr_cal_q,partial_f,dark_frac,ice_clim_f,wse_r_u,slope_r_u,reach_q_b,p_low_slp,p_length,xtrk_dist,obs_frac_n,p_n_nodes')
  if (nzchar(api_key)) {
    # do something
    response = GET(website, add_headers("x-hydrocon-key" = api_key))
  } else {
    # do something else
    response = GET(website)
  } 
  pull = content(response, as='parsed')$results
  if (!is.null(pull)) {
    data = try(read.csv(textConnection(pull$csv), sep=','))
  } else {
    return(NA)
  }
  
  if(is.error(data)){return(NA)}

  data$reach_id = feature_id
  return(data)
}

batch_download_SWOT <- function(obs_ids, api_key){
  print("batch downloading swot")
  # Batch download swot river reaches 
  plan(multisession, workers = workers_)
  SWOT_data = future_lapply(unique(obs_ids), pull_data, api_key = api_key)
  plan(sequential)
  return(SWOT_data)
}

tukey_test = function(ts){
  # print("running tukey")
  # Turkey test for removing outliers
  wseIQR = quantile(ts$wse, c(.25, .75))
  wseT_l = wseIQR[1] - (diff(wseIQR)*1.5)
  wseT_u = wseIQR[2] + (diff(wseIQR)*1.5)
  
  slopeIQR = quantile(ts$slope2, c(.25, .75))
  slopeT_l = slopeIQR[1] - (diff(slopeIQR)*1.5)
  slopeT_u = slopeIQR[2] + (diff(slopeIQR)*1.5)
  
  widthIQR = quantile(ts$width, c(.25, .75))
  widthT_l = widthIQR[1] - (diff(widthIQR)*1.5)
  widthT_u = widthIQR[2] + (diff(widthIQR)*1.5)
  
  ts_filt = ts[ts$width>=widthT_l&ts$width<=widthT_u&ts$wse>=wseT_l&ts$wse<=wseT_u&ts$slope2>=slopeT_l&ts$slope2<=slopeT_u,]
  return(ts_filt)
}

tukey_test_lake = function(ts){
  # print("running tukey lake")
  # Turkey test for removing outliers 
  wseIQR = quantile(ts$wse, c(.25, .75))
  wseT_l = wseIQR[1] - (diff(wseIQR)*1.5)
  wseT_u = wseIQR[2] + (diff(wseIQR)*1.5)
  
  ts_filt = ts[ts$wse>=wseT_l&ts$wse<=wseT_u,]
  return(ts_filt)
}

filter_function = function(swot_ts){
  # print("filtering")
  # Function to filter swot reach data
  # Allowing partial obs to see if it degrades performances. 
  #dawg_filter = tukey_test(swot_ts[swot_ts$time_str!='no_data'&swot_ts$ice_clim_f<2&swot_ts$dark_frac<=0.5&swot_ts$xovr_cal_q<2&!is.na(swot_ts$slope2)&!is.infinite(swot_ts$slope2)&swot_ts$slope2>0&swot_ts$width>0,])
  dawg_filter = tukey_test(swot_ts[swot_ts$time_str!='no_data'&swot_ts$ice_clim_f<2&swot_ts$dark_frac<=0.5&swot_ts$xovr_cal_q<2&swot_ts$partial_f==0&!is.na(swot_ts$slope2)&!is.infinite(swot_ts$slope2)&swot_ts$slope2>0&swot_ts$width>0,])
  perm_relax_filter = tukey_test(swot_ts[swot_ts$time_str!='no_data' & 
                                           swot_ts$p_width >= 60 & 
                                           swot_ts$p_low_slp <= 0 & 
                                           swot_ts$p_length >= 5000 & 
                                           ((swot_ts$xtrk_dist >= 10000 & swot_ts$xtrk_dist <= 60000) | (swot_ts$xtrk_dist >= -60000 & swot_ts$xtrk_dist <= -10000)) & 
                                           swot_ts$ice_clim_f <= 1 & 
                                           swot_ts$reach_q <=2 & 
                                           swot_ts$dark_frac <= 0.6 & 
                                           swot_ts$obs_frac_n > 0.4 & 
                                           swot_ts$xovr_cal_q <= 1 & 
                                           swot_ts$p_n_nodes >= 10,])
  qual_filter = swot_ts[swot_ts$time_str!='no_data'&swot_ts$reach_q<=2,]
  ssf_filter = tukey_test(swot_ts[swot_ts$time_str!='no_data'&swot_ts$reach_q_b<=32768&swot_ts$dark_frac<=0.1&swot_ts$wse_r_u<=0.5&swot_ts$slope_r_u<=10e-5&swot_ts$ice_clim_f==0&swot_ts$xovr_cal_q<=1,])
  return(perm_relax_filter)
}

# Ryan's updated code to get data from farther up/downstream reaches
combining_lk_rv_obs = function(lake, api_key){
  # print("combining lake and river obs")

  #Pull in SWOT river data and subset predownloaded SWOT lake data. 
  upID = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id==lake], ','))
  dnID = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id==lake], ','))
  upObs_all = swot_river[swot_river$reach_id%in%upID,]
  dnObs_all = swot_river[swot_river$reach_id%in%dnID,]
  dnObs_all$reach_id <- as.character(dnObs_all$reach_id)

  sword_continents = list("af", "eu", "as", "as", "oc", "sa", "na", "na", "na")

  #Allow downstream reaches to shift one reach downstream.
  n_ds_reaches = updated_pld$D_reach_n[updated_pld$lake_id==lake]
  n_ds_reaches_obs = dnObs_all[,.N,by=reach_id]
  missing_dn = dnID[dnID%!in%n_ds_reaches_obs$reach_id]
  shift_a_reach_away = function(f, api_key){
    print("Shift a reach away...")

    reach_id = as.character(f)
    reach_continent = strtoi(substring(reach_id, 1, 1))
    continent = paste0(sword_continents[reach_continent], "_sword_v", SWORD_VERSION, ".nc")
    
    sword_nc = nc_open(file.path(indir, "sword", continent))
    sword = list(reach_id=ncvar_get(sword_nc, "reaches/reach_id"),
                 rch_id_dn=ncvar_get(sword_nc, "reaches/rch_id_dn"),
                 n_rch_dn=ncvar_get(sword_nc, "reaches/n_rch_down")
    )
    nc_close(sword_nc)

    n_alt_ds_reaches = sword$n_rch_dn[sword$reach_id==f]

    if(n_alt_ds_reaches==1){
      alt_ds_reach = sword$rch_id_dn[sword$reach_id==f]
      dn_river_pull = try(lapply(alt_ds_reach, pull_data, api_key = api_key))
      dn_river_filt = lapply(dn_river_pull[!is.na(dn_river_pull)], filter_function)
      dnObs_alt = rbindlist(dn_river_filt)
      dnObs_alt$time=as_datetime(dnObs_alt$time_str)
      dnObs_alt$reach_id_original = dnObs_alt$reach_id
      dnObs_alt$reach_id = f
      dnObs_alt$shifted = 'yes'
      if(nrow(dnObs_alt)>0){
        return(dnObs_alt)
      }
    }
  }

  additional_ds_obs = rbindlist(lapply(missing_dn, shift_a_reach_away, api_key = api_key))

  if (nrow(dnObs_all)==0) {
    dnObs_all= additional_ds_obs
  } else {
    additional_ds_obs$reach_id <- as.character(additional_ds_obs$reach_id)
    dnObs_all$reach_id <- as.character(dnObs_all$reach_id)
    dnObs_all = bind_rows(dnObs_all, additional_ds_obs)
  }
  dn_shifted = any(dnObs_all$shifted=='yes')

  #Allow upstream reaches to shift one reach upstream.
  n_us_reaches = updated_pld$D_reach_n[updated_pld$lake_id==lake]
  n_us_reaches_obs = upObs_all[,.N,by=reach_id]
  missing_up = upID[upID%!in%n_us_reaches_obs$reach_id]
  shift_a_reach_up = function(f, api_key){
    print("Shift a reach up")

    reach_id = as.character(f)
    reach_continent = strtoi(substring(reach_id, 1, 1))
    continent = paste0(sword_continents[reach_continent], "_sword_v", SWORD_VERSION, ".nc")

    sword_nc = nc_open(file.path(indir, "sword", continent))
    sword = list(reach_id=ncvar_get(sword_nc, "reaches/reach_id"),
                 rch_id_up=ncvar_get(sword_nc, "reaches/rch_id_up"),
                 n_rch_up=ncvar_get(sword_nc, "reaches/n_rch_up")
                 )
    nc_close(sword_nc)

    n_alt_us_reaches = sword$n_rch_up[sword$reach_id==f]
    if(n_alt_us_reaches==1){
      alt_us_reach = sword$rch_id_up[sword$reach_id==f]
      up_river_pull = try(lapply(alt_us_reach, pull_data, api_key = api_key))
      up_river_filt = lapply(up_river_pull[!is.na(up_river_pull)], filter_function)
      upObs_alt = rbindlist(up_river_filt)
      upObs_alt$time=as_datetime(upObs_alt$time_str)
      upObs_alt$reach_id_original = upObs_alt$reach_id
      upObs_alt$reach_id = f
      upObs_alt$shifted = 'yes'
      if(nrow(upObs_alt)>0){
        return(upObs_alt)
      }
    }
  }

  additional_us_obs = rbindlist(lapply(missing_up, shift_a_reach_up, api_key = api_key))
  additional_us_obs$reach_id <- as.character(additional_us_obs$reach_id)
  upObs_all$reach_id <- as.character(upObs_all$reach_id)

  
  if (nrow(upObs_all)==0) {
    upObs_all= additional_us_obs
  } else {
    upObs_all = bind_rows(upObs_all, additional_us_obs)
  }
  up_shifted = any(upObs_all$shifted=='yes')


  lakeObs_all = lakeFilt[lakeFilt$lake_id==lake,]
  lakeObs_all$time = as_datetime(lakeObs_all$time_str)

  # FIXME: changing lake areas to pld mean lake areas due to SWOT errors. 
  prior_area = updated_pld$Lake_area[updated_pld$lake_id==lake]
  lakeObs_all$area_total = prior_area
  
  if(nrow(lakeObs_all)<3){
    print('noo lake obs... exiting')
    return(NA)}
  if(nrow(upObs_all)<1){
    print('nooo up obs... exiting')
    return(NA)}
  if(nrow(dnObs_all)<1){
    print('nooo down obs... exiting')
    return(NA)}


  ################################################################################
  # Get dates in proper format and subset to matching dates. 
  ################################################################################
  lakeObs_all$date = as.Date(lakeObs_all$time)
  upObs_all$date = as.Date(upObs_all$time)
  dnObs_all$date = as.Date(dnObs_all$time)

  # FIXME: Aggregating lakes to mean values for multiple observations in one day. 
  lakeObs = data.table(lakeObs_all)[,c('wse', 'area_total', 'date')][,lapply(.SD, mean), by=date]
  upObs = data.table(upObs_all)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]
  dnObs = data.table(dnObs_all)[,c('wse', 'width', 'slope', 'slope2','reach_id', 'date')][,lapply(.SD, mean), by=list(date, reach_id)]

  lkDates = unique(lakeObs$date)
  upDts = upObs[,.N,by=date][N>=length(upID)] # limit to dates with obs for each upstream reach.
  dnDts = dnObs[,.N,by=date][N>=length(dnID)] # limit to dates with obs for each downstream reach. 

  #goodDates = lkDates[lkDates%in%upObs_all$date&lkDates%in%dnObs_all$date]
  goodDates = lkDates[lkDates%in%upDts$date&lkDates%in%dnDts$date]

  lakeObsGood = lakeObs[lakeObs$date%in%goodDates,]
  upObsGood = upObs[upObs$date%in%goodDates,]
  dnObsGood = dnObs[dnObs$date%in%goodDates,]

  lakeObs = lakeObsGood[order(lakeObsGood$date),]
  upObs = upObsGood[order(upObsGood$date),]
  dnObs = dnObsGood[order(dnObsGood$date),]

  upObs = upObs[order(upObs$reach_id),]
  dnObs = dnObs[order(dnObs$reach_id),]

  upObs$shifted = up_shifted
  dnObs$shifted = dn_shifted

  if(nrow(lakeObs)<4){return(NA)}
  output = list(lakeObs, upObs, dnObs)
  return(output)
}

download_tributary = function(reaches, start_date='01-01-2023'){
  # print("downloading tribs")
  # Function to pull tributary inflow Q estimates from Geoglows. 
  source_python('src/geoglows_aws_pull.py')
  tributary_flow = pull_tributary(reach_id = reaches, start_date=start_date)
  tributary_flow$date = as.Date(row.names(tributary_flow))
  n_col = ncol(tributary_flow)
  tributary_aggregated = data.table(tributary_flow)[,tributary_total:=rowSums(.SD), .SDcols=-n_col]
  return(tributary_aggregated)
}


pull_geoglows = function(reaches, start_date='01-01-2023'){
  # print("pulling geoglows")
  # Function to pull geoglows data for prior purposes
  source_python('src/geoglows_aws_pull.py')
  model_flow = pull_tributary(reach_id = reaches, start_date=start_date)
  model_flow$Date = as.Date(row.names(model_flow))
  return(data.table(model_flow))
}


extract_data_by_lake <- function(lake, indir){
  # print("extracting data by lake")
    
    # Use dynamic prior Q. False = SOS prior estimate from GRADES / MAF geoglows
    use_ts_prior=TRUE
  
    # Use modeled daily tributary flows. False = mean monthly grades tributaries / MAF geoglows
    use_ts_tributary=TRUE
  
    index=which(names(viable_data)==lake)
    relevant_data = viable_data[index][[1]]
    lakeObs = relevant_data[[1]]
    upObs = relevant_data[[2]]
    dnObs = relevant_data[[3]]

    # Extract month - used to assign et  
    lakeObs$month = month(lakeObs$date)

    # Add the ET data 
    et_lake = et[et$lake_id==lake,]
    if(nrow(et_lake)==0){
        lakeObs$et = 0
    }else{
      et_lake$month = as.numeric(et_lake$month)
      lakeObs$et = et_lake$mean[match(lakeObs$month, et_lake$month)]
    }

    # add in tributary data:Either use geoglow (ts==TRUE or use GRADES-hydroDL mean monthly vals)
    if(use_ts_tributary==TRUE){
        tributary_locations = tributary[tributary$lake_id==(lake),]
        if(nrow(tributary_locations)==0){
            lakeObs$tributary_total=0
        }else{
            tributary_reaches = unique(tributary_locations$LINKNO[tributary_locations$lake_id==lake])
            tributary_reaches = as.list(tributary_reaches)
            tributary_data = download_tributary(tributary_reaches,'01-01-1940')
            mean_annual = mean(tributary_data[,year:=lubridate::year(date)][,mean(tributary_total),year]$V1)
            lakeObs$tributary_total = tributary_data$tributary_total[match(lakeObs$date, tributary_data$date)]
            }
    }else{
        tributary_locations = tributary[tributary$lake_id==(lake),]
        if(nrow(tributary_locations)==0){
            lakeObs$tributary_total=0
        }else{
        tributary_reaches = unique(tributary_locations$LINKNO[tributary_locations$lake_id==lake])
        tributary_reaches = as.list(tributary_reaches)
        tributary_data = download_tributary(tributary_reaches,'01-01-1940')
        mean_annual = mean(tributary_data[,year:=lubridate::year(date)][,mean(tributary_total),year]$V1)
        lakeObs$tributary_total = mean_annual
        }
    }

    # Pull in modeled geoglows data  
    upID = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id==lake], ','))
    dnID = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id==lake], ','))
    sword_reaches = c(upID, dnID)
    sword_geoglows_filt = sword_geoglows[sword_geoglows$reach_id%in%sword_reaches,c('reach_id','LINKNO')]
    geoglows_reaches = unique(as.list(sword_geoglows$LINKNO[sword_geoglows$reach_id%in%sword_reaches]))
    model_data = pull_geoglows(geoglows_reaches, '01-01-1940')

    # Save 
    fwrite(lakeObs, file.path(indir, paste0("clean/lakeobs_", lake, ".csv")))
    fwrite(upObs, file.path(indir, paste0("clean/upobs_", lake, ".csv")))
    fwrite(dnObs,file.path(indir, paste0("clean/dnobs_", lake, ".csv") ))
    fwrite(model_data,file.path(indir,paste0("clean/geoglows_", lake, ".csv")) )

    return()
}




api_key <- get_api_key(prefix)
################################################################################
# Read in lake data via hydrocron
################################################################################

# Load csv of lake ids as a data.table
if (basename(opts$input_file) != "reaches_of_interest.json") {
  
  lakes_input <- fread(opts$input_file)
  lakes_input$lake <- as.character(lakes_input$lake_id)
  
}else{
  reaches_of_interest <- fromJSON(opts$input_file)
  # lake_ids_subset <- lake_ids[grepl("3$", lake_ids)]
  print("reaches of interest")
  print(reaches_of_interest)
  connected_lake_ids <- get_connected_lake_ids(reach_ids = reaches_of_interest, pld = updated_pld)
  lakes_input <- data.table(lake = connected_lake_ids)
  print(lakes_input)
}

files_filt = batch_download_SWOT_lakes(updated_pld$lake_id[updated_pld$lake_id%in%lakes_input$lake], api_key)
combined = rbindlist(files_filt[!is.na(files_filt)])

################################################################################
# Filter lake data
################################################################################
# Testing to see if partial flags make a difference since we're only using wse for lakes at the moment. - partial wse seems great. 
#lakeData = combined[combined$ice_clim_f<2&combined$dark_frac<=0.5&combined$xovr_cal_q<2&combined$time_str!='no_data'&combined$wse>(5000*-1),]
lakeData = combined[combined$ice_clim_f <= 1 &
                      combined$dark_frac <= 0.6 &
                      combined$xovr_cal_q <= 1 &
                      combined$quality_f <= 2 &
                      combined$partial_f <= 0 &
                      ((combined$xtrk_dist >= 10000 & combined$xtrk_dist <= 60000) | (combined$xtrk_dist >= -60000 & combined$xtrk_dist <= -10000)) &
                      combined$time_str!='no_data' &
                      combined$wse>(5000*-1),]
lakeData = lakeData%>%distinct(.keep_all=TRUE)
lakeData$lake_id = as.character(lakeData$lake_id)

# Remove lakedata with multiple ids
lakeData$lake_id_first = sub(";.*", "", lakeData$lake_id)
lakeData$lake_id=lakeData$lake_id_first
lakeData = lakeData[lakeData$lake_id%in%updated_pld$lake_id,]
rm(combined)

################################################################################
# Pull in relevant SWOT river reach data and filter using above functions 
################################################################################
# Filter to lakes with at least n observations. 
n = 5
lakes = unique(data.table(lakeData)[,.N,by=lake_id][N>=n]$lake_id)
lakeFilt = lakeData[lake_id%in%lakes,tukey_test_lake(.SD),by=lake_id]
lakes = unique(data.table(lakeFilt)[,.N,by=lake_id][N>=n]$lake_id)
lakeFilt = lakeFilt[lake_id%in%lakes,]
up_reaches = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id%in%lakes], ','))
dn_reaches = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id%in%lakes], ','))
reaches = c(up_reaches, dn_reaches)
if (is.null(reaches) || length(reaches) == 0){
  message("No reaches found, exiting....")
  quit(save = "no", status = 0)
}
swot_river_pull = batch_download_SWOT(reaches, api_key)
swot_river_filt = lapply(swot_river_pull[!is.na(swot_river_pull)], filter_function)
swot_river = rbindlist(swot_river_filt)
swot_river$time=as_datetime(swot_river$time_str)

# subset to lakes with enough lake and river SWOT obs to run LakeFlow. 
viable_data = lapply(lakes, combining_lk_rv_obs, api_key = api_key)
names(viable_data) = lakes
n_obs_lake = data.table(lake=lakes,obs=unlist(lapply(viable_data, length)))
viable_locations = n_obs_lake[obs>=3,] #Note this 3 is for ensuring an inflow, lake, and outlfow have viable data. 

################################################################################
# Go lake by lake and attach ET and tributary data and save datasets
################################################################################
if (nrow(viable_locations) == 0){
  message("Did not find viable lake, exiting....")
  quit(save = "no", status = 0)
}

for(i in 1:nrow(viable_locations)){
  # print(viable_locations)
    extract_data_by_lake(viable_locations$lake[i], indir)
}

################################################################################
# Save a list of the lakes with data where we will run LakeFlow
################################################################################
# Write out a list of the viable locations
dir.create(file.path(indir, "viable"), showWarnings = FALSE)
numbers <- gregexpr("[0-9]+", basename(opts$input_file))
result <- unlist(regmatches(basename(opts$input_file), numbers))
fwrite(viable_locations[,"lake"], file.path(indir, paste0("viable/viable_locations", index, ".csv")))
print('Found viable lakes...')
