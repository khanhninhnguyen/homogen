#' Three functions in data segmentation.
#' 
#' convert_txt: read raw file daily IWV, used for the segment function
#' forme of raw file : 
#' yymmdd.HHMMSS,  delta_IWV(GPS-ERA),  IWV_GPS,  sigma_ERA5,   sigma_GPS
#' 
#' segment: 
#' 1 total output file 
#' @param path_txt path to raw data in txt 
#' @param list_file list of all files that will be segmented 
#' @param path_result where to store results
#' @param saved_seg_ind do you want to save result of each station?
#' 
#' extract from segmentation results:
#' 2 output files: 
#' * dataname.date_mean.txt: 
#'  station name,   date of begin,   date of end,   mean of this segment 
#' * dataname.monthly_var_stdf.txt
#'  station name,   variance of each month,   std of functional  
#' * dataname.coeff.txt
#' station name,   4 Fourier coefficients 
#' 

library("GNSSfast")
library("gfpop")

#' return NA if it is not txt file 
convert_txt <- function(file_path, saved_folder, saved = 0) {
  data_raw <- read.table(
    file = file_path,
    sep = "",
    skip = 1,
    header = FALSE,
    na.strings = NaN,
    colClasses = c("character", "numeric", "numeric", "numeric", "numeric"),
    col.names = c("date", "delta_IWV", "IWV_GPS", "sigma_ERA5", "sigma_GPS")
  )
  
  name_station <- unlist(strsplit(file_path, "\\."))
  if (name_station[length(name_station)] != "txt") {
    print("Le fichier n'est pas un .txt")
    return(NA)
  }
  
  name_station <- unlist(strsplit(file_path, "\\/"))
  name_station <- substr(name_station[length(name_station)], 13, 16)
  
  ind <- which(substr(data_raw$date, 1, 1) == "9")
  if (length(ind) == 0) {
    data_raw$date <- paste0("20", data_raw$date)
  } else {
    data_raw$date[ind] <- paste0("19", data_raw$date[ind])
    data_raw$date[-ind] <- paste0("20", data_raw$date[-ind])
  }
  
  date_num <- strptime(data_raw$date, "%Y%m%d.%H", tz = "GMT")
  class(date_num)
  unclass(head(date_num, 2))
  
  Y <- data.frame(
    name_station = name_station,
    date = date_num,
    signal = data_raw$delta_IWV,
    ERAI = data_raw$IWV_GPS - data_raw$delta_IWV,
    GPS = data_raw$IWV_GPS
  )
  Y$month <- as.factor(format(Y$date, format = "%m"))
  Y$year <- as.factor(format(Y$date, format = "%Y"))
  
  if(saved != 0){
    file_name = paste0(saved_folder,"/",name_station,".RData")
    save(Y, file = file_name)
  }
  
  return(Y)
}

segment <- function(path_txt, list_file, path_result, criterion, saved_RData = 0, saved_seg_ind = 0, ...){
  
  seg_results <- list()
  extr_info <- list()
  for (i in 1:length(list_file)) {
    
    name_station <- unlist(strsplit(list_file[i], "\\/"))
    name_station <- substr(name_station[length(name_station)], 13, 16)
    
    data.trans <- convert_txt(file_path = paste0(path_txt, list_file[i]),
                              saved_folder = path_result,
                              saved = saved_RData)
    
    GNSS <- GNSSfast(data.trans, lyear = 365.25, lmin = 1, Kmax = 30, 
                     selection.K = "All", S = 0.75, f = TRUE, 
                     selection.f = FALSE, threshold = 0.001, tol = 1e-4)
    
    Tmu = data.frame(begin = data.trans$date[GNSS$seg[[criterion]]$begin], 
                     end = data.trans$date[GNSS$seg[[criterion]]$end], 
                     mean = GNSS$seg[[criterion]]$mean)
    
    out_put <- list( Tmu = Tmu,
                     Stdf = GNSS$funct[[criterion]],
                     Var = GNSS$variances,
                     Coeff = GNSS$coeff[[criterion]])
    
    # save the total results for further studies 
    GNSS$funct <- NULL
    GNSS$Tot <- lapply(GNSS$Tot, function(sublist) {
      sublist[!(names(sublist) %in% "f")]
    })
    
    if(saved_seg_ind != 0){
      save(GNSS, file = paste0(path_result, "tot_seg_", name_station , ".RData"))
      save(out_put, file = paste0(path_result, "seg_", name_station , ".RData"))
    }
    seg_results[[name_station]] <- GNSS
    extr_info[[name_station]] <- out_put
  }
  
  save(seg_results, file = paste0(path_result, "segmentation_total.RData"))
  save(extr_info, file = paste0(path_result, "segmentation_extracted.RData"))
  
}

### extract specific info in txt file   CONTINUE TOMORROW
extract_txt <- function(file.result){
  extract_res = extr_info
  date_mean <- data.frame()
  monthly_var_stdf <- data.frame()
  fourier_coef <- data.frame()
  
}


