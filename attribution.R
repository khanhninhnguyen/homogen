#' 5 functions in change-points attribution.
#' INPUT:
#' * 6 series of difference, 
#' * segmentation results of main and nearby station
#'  (date_mean file from segmentation.R) 
#' 
path_data = "/home/knguyen/data/Data/NGL-attribution/NGL/"
screening <- function(path_data, path_seg_res, path_result, 
                      window_length = 60, min_nb_points, max_coincide){
  
  list_all_pairs = unique(substr(list.files(path_data),6,14))
  
  list_brp = extract_list_brp(a)
  list_main = unique(substr(list_all_pairs,1,4))
  
  # 
}