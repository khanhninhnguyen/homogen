#' 5 functions in change-points attribution.
#' INPUT:
#' * 6 series of difference, 
#' * segmentation results of main and nearby station
#'  (date_mean file from segmentation.R) 
#' 
path_data = "/home/knguyen/data/Data/NGL-attribution/NGL/"
list_six_diff = c("GPS-ERA", "GPS-GPS'", "GPS-ERA'", "ERA-ERA'", "GPS'-ERA'", "GPS'-ERA")
name_six_diff = c("GPS_ERA", "GPS_GPS1", "GPS_ERA1", "ERA_ERA1", "GPS1_ERA1", "GPS1_ERA")

# order: GPS-ERA, GPS-GPS', GPS-ERA', ERA-ERA', GPS'-ERA', GPS'-ERA
  
# DO THIS LATTER 
# screening <- function(path_data, path_seg_res, path_result, 
#                       window_length = 60, min_nb_points, max_coincide){
#   
#   list_all_pairs = unique(substr(list.files(path_data),6,14))
#   
#   list_brp = extract_list_brp(a)
#   list_main = unique(substr(list_all_pairs,1,4))
#   
# 
#   # 
# }


# SELECT NEARBY STATIONS --------------------------------------------------
#' 3 tests: 
#' * limit by other brps 
#' * keep very close changepoints
#' * gap percentage
pre_select_nb <- function(path_data, seg_result, max_limit, min_limit, gap_per){
  
  list_all_pairs = unique(substr(list.files(path_data),6,14))
  list_main = substr(list_all_pairs,1,4)
  list_nearby = substr(list_all_pairs,6,9)
  
  list_brp = extract_list_brp(a)
  
  out = data.frame()
  
  main_st = "rock"
  nearby_sts <- list_nearby[which(list_main == main_st)]
  main_brp = list_brp$brp[which(list_brp$name == main_st)][1]
  main_brps = list_brp$brp[which(list_brp$name == main_st)]
  
  nearby_st = nearby_sts[1]
  nearby_brps = list_brp$brp[which(list_brp$name == nearby_st)]
  
  # merge data 
  df_data = read_data_new(path_data = path_data,
                            main_st = main_st, 
                            nearby_st = nearby_st,
                            name_six_diff = name_six_diff)
  df_data$Date <- as.Date(df_data$Date, format = "%Y-%m-%d")
  # Test 1 

  test_1 = test1(main_brp, nearby_brps)
  
  # Test 2 
  
  test_2 = test2(main_brp, main_brps, df_data, nearby_brps, test_1 = test_1)
  
  test_3 = test3(df_data = df_data, test_1 = test_1, 
                 main_beg = test_2$main_beg_new, 
                 main_end = test_2$main_end_new,
                 nearby_beg = test_2$nearby_beg_new,
                 nearby_end = test_2$nearby_end_new)
}



# CHARACTERIZATION --------------------------------------------------------


