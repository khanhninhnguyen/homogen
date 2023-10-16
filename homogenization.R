rm(list = ls())
# Main steps in homogenization 

## LIB
library("GNSSfast")

## set path: modify the main path, data path and create:
## paths to: store code, result, metadata (if needed) and validation 

path_main = "/home/knguyen/Documents/PhD/"
path_data = "/home/knguyen/data/Data/"
path_results = paste0(path_main,"Results/")
path_code = paste0(path_main,"Code/homogen/")
path_meta = paste0(path_main,"Data/support/")
path_validation = paste0(path_results, "validation/")

# SEGMENTATION ------------------------------------------------------------

## modify to choose list of stations

files = list.files(path_NGL_ERA5)
name_full = substr(files,start = 1, stop = 4)

## need to convert from txt to RData



# VALIDATION --------------------------------------------------------------


# ATTRIBUTION -------------------------------------------------------------


# CORRECTION --------------------------------------------------------------


