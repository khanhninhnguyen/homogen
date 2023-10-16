# Segment data
library("GNSSfast")


# Convert txt to Rdata ----------------------------------------------------
rm(list = ls())

path_main="/home/knguyen/Documents/PhD/"
path_code=paste0(path_main,"Code/convert_data_into_r")
path_data=paste0(path_main,"Data/IWV-ERA/COMPARE_ERA5/")
path_raw=paste0(path_main,"Data/IWV-ERA/NGL_6048/data/")
path_results=paste0(path_main,"Results/")
setwd(path_main)
source(paste0(path_code,"/conversiontoRdata.R"))
files=list.files(paste0(path_raw))
Station.name=substr(files[],1,nchar(files[])-4)

for(i in 1:length(files)){
  
  file = files[i]
  print(file)
  file_name = paste0(path_raw,file)
  #  print(file_name)
  name_folder = "NGL_ERA5.R"
  extract_ngl(file_name, name_folder, path_file = path_data)
  
}

# Run segmentation --------------------------------------------------------

homo_list <- function(station, k, path_data, path_save_homo, nb_test) {
  # k=0: mean running segmentation with 4 criteria, k!=0: no criteria
  for (i in 1:length(station)) {
    print(i)
    station_name <- station[i]
    Y <- get(load(paste0(path_data, station_name, ".RData")))
    if (k == 0) {
      GNSS <- GNSSfast(Y, lyear = 365.25, lmin = 1, Kmax = 30, 
                       selection.K = "All", S = 0.75, f = TRUE, 
                       selection.f = FALSE, threshold = 0.001, tol = 1e-4)
    } else {
      GNSS <- GNSSfast(Y, lyear = 365.25, lmin = 1, Kmax = k, 
                       selection.K = "none", S = 0.75, f = TRUE,
                       selection.f = FALSE, threshold = 0.001, tol = 1e-4)
    }
    save(GNSS, file = paste0(path_save_homo, "homo_", station_name, nb_test, ".RData"))
  }
}



