rm(list = ls())
# Main steps in homogenization 

## LIB
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyverse)
library(nlme)

### uncomment to run segmentation 

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
# ind = c(1:100) # modify to choose list of station to be segmented 
path_raw = paste0("/home/knguyen/data/Data/NGL/data/") 
file_name = "scatterIWV3_0abi.txt"
# file_path = paste0(path_NGL_ERA5, file_name)

segment(path_txt = path_raw,
        criterion = "BM_BJ",
        list_file = list.files(path_raw)[1:2],
        path_result = path_results)


# VALIDATION --------------------------------------------------------------
files = list.files(path_NGL_ERA5)
name_full = substr(files,start = 1, stop = 4)

# ATTRIBUTION -------------------------------------------------------------


# CORRECTION --------------------------------------------------------------


