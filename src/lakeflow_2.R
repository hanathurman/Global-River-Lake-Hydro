################################################################################
# Run Lakeflow
################################################################################

################################################################################
# Load modules 
################################################################################
library(RNetCDF)
library(foreign)
library(lubridate)
library(rstan)
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
'%!in%' <- function(x,y)!('%in%'(x,y))

# Example Deployment

# docker run -v /mnt/lakeflow:/data/input -v /mnt/input/sos:/sos/input lakeflow_deploy -c /data/input/viable/viable_locations1.csv -s /sos/input -w 1 -i /data/input -o /data/input/out

################################################################################
# Set args
option_list <- list(
    make_option(c("-c", "--input_file"), type = "character", default = NULL, help = "filepath for a list of lakes to run lakeflow"),
    make_option(c("-s", "--sos_dir"), type = "character", default = NULL, help = "filepath for the sos"),
    make_option(c("-w", "--workers"), type = "integer", default = NULL, help = "number of cores for lakeflow to use"),
    make_option(c("-i", "--indir"), type = "character", default = NULL , help = "directory with input files"),
    make_option(c("-o", "--outdir"), type = "character", default = NULL , help = "directory to output results"),
    make_option(c("--index"), type = "integer", default = NULL, help = "index of what lake to process in the input file if blank process whole file")
  )
################################################################################

# Grab arguments
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# Load a df of lakes where lakeflow will be run
lakes <- fread(opts$input_file)
#lakes <- fread("viable_locations.csv")
lakes$lake <- as.character(lakes$lake)

#lakes <- data.table(lake=c("1140003043"))

# Set the number of cores that lakeflow will use
cores = as.numeric(opts$workers)
#cores <- 6

# Path that holds input data, was previously 'in'
indir <- opts$indir

# Path that holds sos
sos_dir <- opts$sos_dir

# Path that you write to
outdir <- opts$outdir

# Set index, use aws array if index is -256 (ascii code for AWS)
index <- opts$index

if (!(is.null(index))){
    if (index == -256){
    index <- strtoi(Sys.getenv("AWS_BATCH_JOB_ARRAY_INDEX")) + 1
    }
    else{
        index <- index + 1
    }
}

################################################################################
# Load datasets 
################################################################################

updated_pld = fread(file.path(indir,"SWORDv16_PLDv103_wo_ghost_rch.csv"))
updated_pld$lake_id =  as.character(updated_pld$lake_id)
updated_pld$continent = substr(updated_pld$lake_id, 1,1)

sword_geoglows = fread(file.path(indir,'ancillary/sword_geoglows.csv'))
sword_geoglows$reach_id = as.character(sword_geoglows$reach_id)

################################################################################
# Define functions 
################################################################################

lakeFlow = function(lake){
    
    # Use dynamic prior Q. False = SOS prior estimate from GRADES / MAF geoglows
    use_ts_prior=TRUE
    
    # Use modeled daily tributary flows. False = mean monthly grades tributaries / MAF geoglows
    use_ts_tributary=TRUE

    lakeObs <- fread(file.path(indir, paste0("clean/lakeobs_", lake, ".csv")))
    upObs <- fread(file.path(indir, paste0("clean/upobs_", lake, ".csv")))
    upObs$reach_id <- as.character(upObs$reach_id)
    dnObs <- fread(file.path(indir, paste0("clean/dnobs_", lake, ".csv")))
    dnObs$reach_id <- as.character(dnObs$reach_id)
    
    upID = unlist(strsplit(updated_pld$U_reach_id[updated_pld$lake_id==lake], ','))
    dnID = unlist(strsplit(updated_pld$D_reach_id[updated_pld$lake_id==lake], ','))
    
    d_x_area = function(wse, width){
        wse_from_min = wse - min(wse)
        x_area = width * wse_from_min
        area_dt = c(NA, diff(x_area))
        return(x_area)
    }
  
    # Apply d_x_area function from above. 
    upObs$d_x_area = upObs[,d_x_area(wse, width),by=reach_id]$V1
    dnObs$d_x_area = dnObs[,d_x_area(wse, width), by=reach_id]$V1
  
    # Calculate lake storage and storage change
    lake_wse_from_min = lakeObs$wse - min(lakeObs$wse)
    area_m2 = lakeObs$area_total * 1e6 # convert from km2 to m2
    lakeObs$storage = area_m2 * lake_wse_from_min
    ht_change = c(NA, diff(lakeObs$wse))
    area_m2 = lakeObs$area_total*1e6

    area_val = list()
    for(k in 2:length(area_m2)){
        current = area_m2[k]
        prior = area_m2[k-1]
        area_add = current+prior
        area_mult= current*prior
        area_sqrt = sqrt(area_mult)
        area_val[[k]] = area_add+area_sqrt
    }
    area_param = c(NA, unlist(area_val))
    lakeObs$storage_dt = (ht_change*area_param)/3
  
    # Put data into table matching LakeFlow code (synthetic dataset): 
    lakeObsOut = lakeObs
    upObsOut = upObs
    dnObsOut = dnObs
    names(lakeObsOut) = c("date_l", "wse_l", "area_total_l", "month", "et", "tributary_total", "storage_l", "storage_dt_l")
    names(upObsOut) = paste0(names(upObs), "_u")
    names(dnObsOut) = paste0(names(dnObs), "_d")
  
    # Split dataframes by group into lists of dataframes.  
    up_df_list = split(upObsOut, by='reach_id_u')
    dn_df_list = split(dnObsOut, by='reach_id_d')
  
    # n days between observations. 
    lakeObsOut$n_days = c(NA, as.numeric(diff(lakeObsOut$date_l)))
  
    # Remove first day from all dataframes bc dv is missing. 
    lakeObsOut = lakeObsOut[-1,]
    up_df_list = lapply(up_df_list, function(f) f[-1,])
    dn_df_list = lapply(dn_df_list, function(f) f[-1,])
  
    # Transpose the data to get in proper format for stan. 
    up_df_wide = lapply(up_df_list, t)
    dn_df_wide = lapply(dn_df_list, t)
  
    # Combine the swot obs of multiple reaches back together.  
    up_df_out = do.call('rbind', up_df_wide)
    dn_df_out = do.call('rbind', dn_df_wide)
  
    # Create seperate matrices for each type of swot obs. Time across, reach down
    up_df_stan = lapply(names(upObsOut)[-1], function(f) as.matrix(up_df_out[grep(f,rownames(up_df_out)),]))
    names(up_df_stan) = names(upObsOut)[-1]
    up_df_stan = lapply(up_df_stan, apply, 2, as.numeric)
  
    dn_df_stan = lapply(names(dnObsOut)[-1], function(f) as.matrix(dn_df_out[grep(f,rownames(dn_df_out)),]))
    names(dn_df_stan) = names(dnObsOut)[-1]
    dn_df_stan = lapply(dn_df_stan, apply, 2, as.numeric)
  
    transposer <- function(df) {
        z<-t(df)
    }
  
    # If reach is only one value, transpose the data to match multiple inflow/outflow setups. 
    if(length(upID)==1){
        up_df_stan = lapply(up_df_stan, transposer)
    }
  
    if(length(dnID)==1){
        dn_df_stan = lapply(dn_df_stan, transposer)
    }
  
    # dA shift function. 
    da_shift_fun = function(x) median(x) - min(x)
  
    # dA scale to be greater than 0.
    da_scale = function(x) x - min(x)
  
    # Apply the above functions. 
    up_da_shift = array(apply(up_df_stan$d_x_area_u,1, da_shift_fun))
    d_x_area_scaled_u = t(as.matrix(apply(up_df_stan$d_x_area_u, 1, da_scale)))
  
    dn_da_shift = array(apply(dn_df_stan$d_x_area_d,1, da_shift_fun))
    d_x_area_scaled_d = t(as.matrix(apply(dn_df_stan$d_x_area_d, 1, da_scale)))
  
    # Function to extract priors from SOS. 
    sos_pull = function(reach_id){
        #sos = "sos/constrained/na_sword_v15_SOS_priors.nc"
        #sos_outflow = RNetCDF::open.nc(sos)
        #sos = "in/sos/constrained/na_sword_v15_SOS_priors.nc"
        if (updated_pld$continent[updated_pld$lake_id==lake] == 1) {
          sos = file.path(sos_dir, "af_sword_v16_SOS_priors.nc")
        } else if (updated_pld$continent[updated_pld$lake_id==lake] == 2) {
          sos = file.path(sos_dir, "eu_sword_v16_SOS_priors.nc")
        } else if (updated_pld$continent[updated_pld$lake_id==lake] == 3) {
          sos = file.path(sos_dir, "as_sword_v16_SOS_priors.nc")
        } else if (updated_pld$continent[updated_pld$lake_id==lake] == 4) {
          sos = file.path(sos_dir, "as_sword_v16_SOS_priors.nc")
        } else if (updated_pld$continent[updated_pld$lake_id==lake] == 5) {
          sos = file.path(sos_dir, "oc_sword_v16_SOS_priors.nc")
        } else if (updated_pld$continent[updated_pld$lake_id==lake] == 6) {
          sos = file.path(sos_dir, "sa_sword_v16_SOS_priors.nc")
        } else {sos = file.path(sos_dir, "na_sword_v16_SOS_priors.nc")} #Sets NA priors for continents 7, 8, and 9
        sos_outflow = RNetCDF::open.nc(sos)
        reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
        reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
        index <- which(reach_ids==reach_id, arr.ind=TRUE)
        gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
        geoA = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
        geoN = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
        geoNlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
        geoNupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
        geoAlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
        geoAupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
        geoAsd = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
        geoNsd = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
        model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
        qHat <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
        #Assign a prior of 1000 cms is q prior is NA. - probably makes sense to use some relationship between width and meanQ in the future.   
        qHat = ifelse(is.na(qHat), 1000, qHat)
        qUpper <- RNetCDF::var.get.nc(model_grp, "max_q")[index]
        qLower <- RNetCDF::var.get.nc(model_grp, "min_q")[index]
        qLower = ifelse(qLower<=0|is.na(qLower), 0.001, qLower)
        qUpper = ifelse(is.na(qUpper),500000, qUpper)
        sigma = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
        qSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
        return(data.table(reach_id, geoA, geoN, geoNlower, geoNupper, geoAlower, geoAupper,
                        geoAsd, geoNsd, qHat, qUpper, qLower, sigma, qSd))
    }

    # Create our own geobam priors when SoS provides null priors - common issue bc lake influenced reaches often don't pass SoS QAQC.
    sos_fit = function(df){
        reach_id = unique(unlist(df[grep('reach_id', names(df))]))#[[1]]
        width_matrix = df[grep('width', names(df))][[1]]
        slope_matrix = df[grep('slope2', names(df))][[1]]
        geoA = estimate_logA0(width_matrix)
        geoN = estimate_logn(slope_matrix, width_matrix)
        geoNlower = estimate_lowerboundlogn(width_matrix)
        geoNupper = estimate_upperboundlogn(width_matrix)
        geoAlower = estimate_lowerboundA0(width_matrix)
        geoAupper = estimate_upperboundA0(width_matrix)
        geoAsd = estimate_A0SD(width_matrix)
        geoNsd = estimate_lognSD(width_matrix)
        #Talk to Craig about more informed estimates of the below. 
        qHat = rep(NA, length(reach_id))
        qLower = rep(NA, length(reach_id))
        qUpper = rep(NA, length(reach_id))
        qSd = rep(100, length(reach_id))
        sigma = rep(0.25, length(reach_id))
        #Assign a prior of 1000 cms is q prior is NA. - probably makes sense to use some relationship between width and meanQ in the future.   
        qHat = ifelse(is.na(qHat), 1000, qHat)
        qLower = ifelse(qLower<=0|is.na(qLower), 0.001, qLower)
        qUpper = ifelse(is.na(qUpper),500000, qUpper)
        return(list(data.table(reach_id, geoA, geoN, geoNlower, geoNupper, geoAlower, geoAupper,
                      geoAsd, geoNsd, qHat, qUpper, qLower, sigma, qSd)))
    }
  
    # Add u and d labels to each column of sos outputs - helps with organizing the data. 
    nms_paste = function(df, added_label){
        nms = colnames(df)
        new_nms = paste0(nms, '_', added_label)
        colnames(df) = new_nms
        return(df)
    }
  
    # Pull priors from sos 
    #New 'any(is.na)...' is new to account for NA values in addition to zeros. 
    up_sos = lapply(upID, sos_pull)
    sos_geobam = 'sos'
    if(any(unlist(lapply(up_sos, nrow))==0)|any(is.na(unlist(up_sos)))){
        up_sos = sos_fit(up_df_stan)
        sos_geobam = 'geobam'
    }
    up_sos = lapply(up_sos, nms_paste, 'u')
  
    dn_sos = lapply(dnID, sos_pull)
    if(any(unlist(lapply(dn_sos, nrow))==0)|any(is.na(unlist(dn_sos)))){
        dn_sos = sos_fit(dn_df_stan)
        sos_geobam = 'geobam'
    }
    dn_sos = lapply(dn_sos, nms_paste, 'd')
  
    # Transpose the sos data to get in proper format for stan. 
    up_sos_wide = lapply(up_sos, t)
    dn_sos_wide = lapply(dn_sos, t)
  
    # Combine the sos priors of multiple reaches together.  
    up_sos_out = do.call('rbind', up_sos_wide)
    dn_sos_out = do.call('rbind', dn_sos_wide)
  
    # Create seperate matrices for each type of sos prior. reach down time across. 
    up_sos_stan = lapply(colnames(up_sos[[1]]), function(f) as.matrix(up_sos_out[grep(f,rownames(up_sos_out)),]))
    names(up_sos_stan) = colnames(up_sos[[1]])
    up_sos_stan = lapply(up_sos_stan, apply, 2, as.numeric)
  
    dn_sos_stan = lapply(colnames(dn_sos[[1]]), function(f) as.matrix(dn_sos_out[grep(f,rownames(dn_sos_out)),]))
    names(dn_sos_stan) = colnames(dn_sos[[1]])
    dn_sos_stan = lapply(dn_sos_stan, apply, 2, as.numeric)
  
    # Place in array to meet stan needs. 
    up_sos_stan = lapply(up_sos_stan, array)
    dn_sos_stan = lapply(dn_sos_stan, array)
  
    # Create matrix for select priors that need to be length of obs for stan. 
    convert_to_matrix = function(f){
        n = length(f)
        output = apply(f, 1, rep, nrow(lakeObsOut))
        return(t(output))
    }
  
    up_sos_stan$sigma_u=convert_to_matrix(up_sos_stan$sigma_u)
    dn_sos_stan$sigma_d=convert_to_matrix(dn_sos_stan$sigma_d)
    up_sos_stan$qSd_u = convert_to_matrix(up_sos_stan$qSd_u)
    dn_sos_stan$qSd_d = convert_to_matrix(dn_sos_stan$qSd_d)
    up_sos_stan$qHat_u = convert_to_matrix(up_sos_stan$qHat_u)
    dn_sos_stan$qHat_d = convert_to_matrix(dn_sos_stan$qHat_d)
  
    if(length(up_sos_stan$reach_id_u)==0){return(NA)}
    if(length(dn_sos_stan$reach_id_d)==0){return(NA)}
  
  
    # Pull in Modeled timeseries as prior Q rather than mean: #Updating alternative to using static Geoglows mean annual flow. 
    if(use_ts_prior==TRUE){
        #sword_geoglows = fread(paste0(inPath, '/ancillary/sword_geoglows.csv'))
        sword_reaches = c(upID, dnID)
        #sword_geoglows$reach_id = as.character(sword_geoglows$reach_id)
        sword_geoglows_filt = sword_geoglows[sword_geoglows$reach_id%in%sword_reaches,c('reach_id','LINKNO')]
        geoglows_reaches = unique(as.list(sword_geoglows$LINKNO[sword_geoglows$reach_id%in%sword_reaches]))
        # Pull in modeled geoglows data. 
        #model_data = pull_geoglows(geoglows_reaches, '01-01-1940')
        model_data <- fread(file.path(indir,paste0("clean/geoglows_", lake, ".csv")))
        model_data$Date <- as.Date(model_data$Date)
        # Convert bad Q to NA. 
        ind = ncol(model_data)-1
        model_data[,1:ind][model_data[,1:ind] < 0] <- 0.1
        # When model data is unavailable (NRT SWOT data), use mean model. 
        missing_dates = data.table(Date=seq.Date((max(model_data$Date)+1), Sys.Date(), by=1))
        model_data = bind_rows(model_data, missing_dates)
        model_data[] <- lapply(model_data, function(x) { 
            x[is.na(x)] <- mean(x, na.rm = TRUE)
            x})
        model_data_all = model_data
        model_data = model_data[model_data$Date%in%lakeObsOut$date_l,]
        model_wide = melt(model_data, id.vars=c('Date'))
        model_wide$LINKNO = as.character(model_wide$variable)
        model_wide$reach_id = sword_geoglows_filt$reach_id[match(model_wide$LINKNO, sword_geoglows_filt$LINKNO)]
        model_list = split(model_wide, by='reach_id')
    
        up_ind = lapply(up_sos_stan$reach_id_u,function(x){which(as.character(x)==names(model_list))})
        dn_ind = lapply(dn_sos_stan$reach_id_d,function(x){which(as.character(x)==names(model_list))})
    
        if(length(unique(model_wide$LINKNO))==1){
            up_ind[[1]] = 1
            dn_ind[[1]] = 1
            }
    
        if(length(up_ind[[1]])==0){
            up_ind[[1]] = 2
            }
    
        up_qhat = lapply(unlist(up_ind), function(x){t(as.matrix(model_list[[x]]$value))})
        up_qhat = do.call(rbind, up_qhat)
    
        dn_qhat = lapply(unlist(dn_ind), function(x){t(as.matrix(model_list[[x]]$value))})
        dn_qhat = do.call(rbind, dn_qhat)
    
        up_sos_stan$qHat_u = up_qhat
        dn_sos_stan$qHat_d = dn_qhat
    
    }else{
        #sword_geoglows = fread(paste0(inPath, '/ancillary/sword_geoglows.csv'))
        sword_reaches = c(upID, dnID)
        #sword_geoglows$reach_id = as.character(sword_geoglows$reach_id)
        sword_geoglows_filt = sword_geoglows[sword_geoglows$reach_id%in%sword_reaches,c('reach_id','LINKNO')]
        geoglows_reaches = unique(as.list(sword_geoglows$LINKNO[sword_geoglows$reach_id%in%sword_reaches]))
        # Pull in modeled geoglows data. 
        #model_data = pull_geoglows(geoglows_reaches, '01-01-1940')
        model_data <- fread(file.path(indir,paste0("clean/geoglows_", lake, ".csv")))
        model_data$Date <- as.Date(model_data$Date)
        # Convert bad Q to NA. 
        ind = ncol(model_data)-1
        model_data[,1:ind][model_data[,1:ind] < 0] <- 0.1
        # When model data is unavailable (NRT SWOT data), use mean model. 
        missing_dates = data.table(Date=seq.Date((max(model_data$Date)+1), Sys.Date(), by=1))
        model_data = bind_rows(model_data, missing_dates)
        model_data[] <- lapply(model_data, function(x) { 
            x[is.na(x)] <- mean(x, na.rm = TRUE)
            x})
        #model_data = model_data[model_data$Date%in%lakeObsOut$date_l,]
        model_annual = model_data[,year:=lubridate::year(Date),]
        model_wide = melt(model_data, id.vars=c('year'))
        model_wide = model_wide[,list(value=mean(value)),by=list(year, variable)][,list(value=mean(value)),by=variable][variable!='Date']
        model_wide = model_wide[rep(model_wide[,.I], nrow(lakeObsOut))]
    
        model_wide$LINKNO = as.character(model_wide$variable)
        model_wide$reach_id = sword_geoglows_filt$reach_id[match(model_wide$LINKNO, sword_geoglows_filt$LINKNO)]
        model_list = split(model_wide, by='reach_id')
    
        up_ind = lapply(up_sos_stan$reach_id_u,function(x){which(as.character(x)==names(model_list))})
        dn_ind = lapply(dn_sos_stan$reach_id_d,function(x){which(as.character(x)==names(model_list))})
    
        if(length(unique(model_wide$LINKNO))==1){
            up_ind[[1]] = 1
            dn_ind[[1]] = 1
        }
    
        up_qhat = lapply(unlist(up_ind), function(x){t(as.matrix(model_list[[x]]$value))})
        up_qhat = do.call(rbind, up_qhat)
    
        dn_qhat = lapply(unlist(dn_ind), function(x){t(as.matrix(model_list[[x]]$value))})
        dn_qhat = do.call(rbind, dn_qhat)
    
        up_sos_stan$qHat_u = up_qhat
        dn_sos_stan$qHat_d = dn_qhat
    }
  
    #New: Stan can't handle 0 prior estimates of discharge. So it replaces it with the mean value. 
    up_sos_stan$qHat_u = ifelse(up_sos_stan$qHat_u<=0, mean(up_sos_stan$qHat_u), up_sos_stan$qHat_u)
    dn_sos_stan$qHat_d = ifelse(dn_sos_stan$qHat_d<=0, mean(dn_sos_stan$qHat_d), dn_sos_stan$qHat_d)
  
    # Combine all data needed for stan. 
    stan_data = list(N = nrow(lakeObsOut),
                    n1=length(upID),
                    n2=length(dnID),
                    sigmaIn = up_sos_stan$sigma_u,
                    sigmaOut = dn_sos_stan$sigma_d,
                    qInSd = up_sos_stan$qSd_u,#rep(qInSd,nrow(data)),#qInSd used to multiply by 0.1?
                    qOutSd = dn_sos_stan$qSd_d,#rep(qOutSd,nrow(data)),#qOutSd used to multiply by 0.1?
                    q = log(up_sos_stan$qHat_u),#rep(log(qHatIn),nrow(data)),
                    sigma = up_sos_stan$sigma_u,#rep(0.25, nrow(data)),
                    da = d_x_area_scaled_u,#data$upArea_dt_u-min(data$upArea_dt_u, na.rm=T), # possible error from Ryan: dA is already change in A, but he's subtracting min dA 
                    w=log(up_df_stan$width_u),
                    s=log(up_df_stan$slope2_u), # using smoothed slope (slope2) to minimize negatives
                    da2=d_x_area_scaled_d,#da2=data$dnArea_dt_d-min(data$dnArea_dt_d, na.rm=T), # possible error from Ryan: dA is already change in A, but he's subtracting min dA 
                    w2=log(dn_df_stan$width_d),
                    s2=log(dn_df_stan$slope2_d), # = smoothed slope (slope2) produced similar results as raw slope 
                    q2=log(dn_sos_stan$qHat_d),#rep(log(qHatOut),nrow(data)),
                    #dv=data$storage_dt_l/86400, #seconds per day
                    et=lakeObsOut$et,#rep(0, nrow(lakeObsOut)), # filled with zero vals right now
                    lateral=lakeObsOut$tributary_total,#rep(0,nrow(lakeObsOut)), # filled with zero vals right now
                    dv_per = (lakeObsOut$storage_dt_l/86400)/lakeObsOut$n_days)
  
  
    # set NAs, NaNs, and Inf to 0 for stan (Ryan's code):
    stan_data$da = ifelse(is.na(stan_data$da)|is.infinite(stan_data$da), 0, stan_data$da)
    stan_data$da2 = ifelse(is.na(stan_data$da2)|is.infinite(stan_data$da2), 0, stan_data$da2)
    stan_data$s = ifelse(is.na(stan_data$s), 0, stan_data$s)
    stan_data$s2 = ifelse(is.na(stan_data$s2), 0, stan_data$s2)
    stan_data$lateral = ifelse(is.na(stan_data$lateral), 0, stan_data$lateral)
    stan_data$et = ifelse(is.na(stan_data$et), 0, stan_data$et)
    stan_data$nInlower = up_sos_stan$geoNlower_u #as.array(geoNInlower)
    stan_data$nInupper = up_sos_stan$geoNupper_u#as.array(geoNInupper)
    stan_data$aInlower= up_sos_stan$geoAlower_u#as.array(geoAInlower)
    stan_data$aInupper = up_sos_stan$geoAupper_u#as.array(geoAInupper)
    stan_data$nOutlower = dn_sos_stan$geoNlower_d#as.array(geoNOutlower)
    stan_data$nOutupper = dn_sos_stan$geoNupper_d#as.array(geoNOutupper)
    stan_data$aOutlower= dn_sos_stan$geoAlower_d#as.array(geoAOutlower)
    stan_data$aOutupper = dn_sos_stan$geoAupper_d#as.array(geoAOutupper)
    stan_data$daInShift = up_da_shift#apply(up_df_stan$d_x_area_u,1, da_shift_fun)#as.array(da_shift_fun(data$upArea_dt_u))
    stan_data$daOutShift = dn_da_shift#apply(dn_df_stan$d_x_area_d,1, da_shift_fun)#as.array(da_shift_fun(data$dnArea_dt_d))
    stan_data$nInHat = up_sos_stan$geoN_u#as.array(geoNin)
    stan_data$nInSd = up_sos_stan$geoNsd_u#as.array(0.25) # Ryan, what's this? 
    stan_data$aInHat = up_sos_stan$geoA_u#as.array(geoAin)
    stan_data$aInSd = up_sos_stan$geoAsd_u#as.array(geoAinSD)
    stan_data$qInupper=log(up_sos_stan$qUpper_u)#as.array(log(qInupper))
    stan_data$qInlower=log(up_sos_stan$qLower_u)#as.array(log(qInlower))
    stan_data$nOutHat = dn_sos_stan$geoN_d#as.array(geoNout) 
    stan_data$nOutSd = dn_sos_stan$geoNsd_d#as.array(0.25)
    stan_data$aOutHat = dn_sos_stan$geoA_d#as.array(geoAout)
    stan_data$aOutSd = dn_sos_stan$geoAsd_d#as.array(geoAoutSD)
    stan_data$qOutupper=log(dn_sos_stan$qUpper_d)#as.array(log(qOutupper))
    stan_data$qOutlower=log(dn_sos_stan$qLower_d)#as.array(log(qOutlower))
  
    # Apply stan code. 
    fit = try(stan("src/lakeflow_stan_flexible.stan",
                    data=stan_data,
                    chains=6,#4, #3
                    cores=cores, #6
                    iter=4000,#4000, #iter=4000
                    control=list(stepsize=0.5,
                                adapt_delta=0.9)))
    if(is.error(fit)){next}
  
    # Manning's eqn:
    eqn1 = function(n, a, da, w, s){
        flow = (n^-1)*((a+da)^(5/3))*(w^(-2/3))*(s^0.5)
        return(flow)
    }
  
    # Pull outputs from stan. 
    mn = get_posterior_mean(fit)
    rw = row.names(mn)
    mn = data.table(mn)
    mn$type = rw
    if(nrow(mn)==0|ncol(mn)<6){return(NA)}
    mn = mn[, c('type','mean-all chains')]
  
    uncertainty = data.table(summary(fit)$summary)
    uncertainty$type = mn$type
  
    # Manning's n estimates from stan
    roughness_inflow = mn[startsWith(mn$type,'n[')]
    roughness_inflow_sd = uncertainty[startsWith(uncertainty$type,'n[')]
    roughness_outflow = mn[startsWith(mn$type,'nOut[')]
    roughness_outflow_sd = uncertainty[startsWith(uncertainty$type,'nOut[')]
  
    # A0 estimates from stan
    bath_inflow = mn[startsWith(mn$type,'a[')]
    bath_inflow_sd = uncertainty[startsWith(uncertainty$type,'a[')]
    bath_outflow = mn[startsWith(mn$type,'aOut[')]
    bath_outflow_sd = uncertainty[startsWith(uncertainty$type,'aOut[')]
  
    # Bayesian estimate of inflow and outflow - not used for our purposes. 
    bayes_inflow = mn[startsWith(mn$type,'logQ_in[')]
    bayes_inflow_sd = uncertainty[startsWith(uncertainty$type,'logQ_in[')]$sd
    bayes_outflow = mn[startsWith(mn$type,'logQ_out[')]
    bayes_outflow_sd = uncertainty[startsWith(uncertainty$type,'logQ_out[')]$sd
  
    # Quick way to organize the data in case of varying N of inflows/outflows - could be placed in a function at some point. 
    inflow_outputs = list()
    for(j in 1:length(upID)){
        q_estimate = eqn1(exp(roughness_inflow$`mean-all chains`[j]),bath_inflow$`mean-all chains`[j],
                      d_x_area_scaled_u[j,],up_df_stan$width_u[j,], up_df_stan$slope2_u[j,])
        q_bayes = bayes_inflow$`mean-all chains`[bayes_inflow$type]
        inflow_outputs[[j]] = data.table(q_lakeflow = q_estimate, reach_id=upID[j], n_lakeflow=exp(roughness_inflow$`mean-all chains`[j]),a0_lakeflow=bath_inflow$`mean-all chains`[j],date=lakeObsOut$date_l,lake_id=lake,q_model=up_sos_stan$qHat_u[j,],
                                     width = stan_data$w[j,], slope2 = stan_data$s[j,], da = d_x_area_scaled_u[j,], wse = up_df_stan$wse_u[j,],
                                     storage=lakeObsOut$storage_l,dv=stan_data$dv_per, tributary=stan_data$lateral, et=stan_data$et,type='inflow',n_lakeflow_sd = roughness_inflow_sd$sd[j],a0_lakeflow_sd =bath_inflow_sd$sd[j],
                                     q_upper = stan_data$qInupper[j],q_lower = stan_data$qInlower[j])
    }
    inflow_outputs = rbindlist(inflow_outputs)
    inflow_outputs$bayes_q = bayes_inflow$`mean-all chains`
    inflow_outputs$bayes_q_sd = bayes_inflow_sd
  
    outflow_outputs = list()
    for(j in 1:length(dnID)){
        q_estimate = eqn1(exp(roughness_outflow$`mean-all chains`[j]),bath_outflow$`mean-all chains`[j],
                      d_x_area_scaled_d[j,],dn_df_stan$width_d[j,], dn_df_stan$slope2_d[j,])
        outflow_outputs[[j]] = data.table(q_lakeflow = q_estimate, reach_id=dnID[j], n_lakeflow=exp(roughness_outflow$`mean-all chains`[j]),a0_lakeflow=bath_outflow$`mean-all chains`[j],date=lakeObsOut$date_l,lake_id=lake,q_model=dn_sos_stan$qHat_d[j,],
                                      width = stan_data$w2[j,], slope2 = stan_data$s2[j,], da = d_x_area_scaled_d[j,], wse = dn_df_stan$wse_d[j,],
                                      storage=lakeObsOut$storage_l,dv=stan_data$dv_per, tributary=stan_data$lateral, et=stan_data$et,type='outflow',n_lakeflow_sd = roughness_outflow_sd$sd[j],a0_lakeflow_sd =bath_outflow_sd$sd[j],
                                      q_upper = stan_data$qOutupper[j],q_lower = stan_data$qOutlower[j])
    }
    outflow_outputs = rbindlist(outflow_outputs)
    outflow_outputs$bayes_q = bayes_outflow$`mean-all chains`
    outflow_outputs$bayes_q_sd = bayes_outflow_sd
  
  
  
    output_df = bind_rows(inflow_outputs, outflow_outputs)
    output_df$prior_fit = sos_geobam
    fwrite(output_df, file.path(outdir, paste0(lake, '.csv')))

    return()

    }


################################################################################
# Apply lakeflow to the list of lakes
################################################################################

if (is.null(index)){
    lake_data <- lakes  
}else{
    lake_data <- lakes[index, , drop = FALSE]
}

for(i in 1:nrow(lake_data)){
  lakeFlow(lake_data$lake[i])
}

