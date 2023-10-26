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
  
  date_mean = read.table(file = paste0("/home/knguyen/Documents/PhD/Results/", "NGL6048.BM_BJ.date_mean.txt"), header = TRUE)
  list_all_pairs = unique(substr(list.files(path_data),6,14))
  list_main = substr(list_all_pairs,1,4)
  list_nearby = substr(list_all_pairs,6,9)
  list_brp = extract_list_brp(date_mean) # after removed cluster of breaks by the first
  list_main_unique = unique(list_main)
  list_main_inhomo = list_main_unique[which(list_main_unique %in% list_brp$name)]
  
  out = data.frame(main = character(0),
                   brp = as.Date(character(0)),     
                   nearby = character(0),
                   coincident = numeric(0),
                   dist_noise = numeric(0),
                   main_beg_new = as.Date(character(0)),      
                   main_end_new = as.Date(character(0)),        
                   nearby_beg_new = as.Date(character(0)),     
                   nearby_end_new = as.Date(character(0)),      
                   n_main_bef = numeric(0),
                   n_nearby_bef = numeric(0),
                   n_joint_bef = numeric(0),
                   n_main_aft = numeric(0),
                   n_nearby_aft = numeric(0), 
                   n_joint_aft = numeric(0),
                   r_main_bef = numeric(0),
                   r_nearby_bef = numeric(0),
                   r_joint_bef = numeric(0),  
                   r_main_aft = numeric(0),
                   r_nearby_aft = numeric(0),   
                   r_joint_aft = numeric(0))

  for (i in c(21:length(list_main_inhomo))) {
    print(i)
    main_st = list_main_inhomo[i]
    main_brps = list_brp$brp[which(list_brp$name == main_st)]
    nearby_sts <- list_nearby[which(list_main == main_st)]
    
    for (j in c(1:length(main_brps))) {
      main_brp = main_brps[j]
      
      for (k in c(1:length(nearby_sts))) {
        nearby_st = nearby_sts[k]
        nearby_brps = list_brp$brp[which(list_brp$name == nearby_st)]
        # merge data 
        df_data = read_data_new(path_data = path_data,
                                main_st = main_st, 
                                nearby_st = nearby_st,
                                name_six_diff = name_six_diff)
        df_data$Date <- as.Date(df_data$Date, format = "%Y-%m-%d")
        # Test 1 
        
        test_1 <- test1(main_brp = main_brp, nearby_brps = nearby_brps)
        
        # Test 2 
        
        test_2 <- test2(main_brp = main_brp,
                       main_brps = main_brps, 
                       df_data = df_data,
                       nearby_brps = nearby_brps, 
                       test_1 = test_1)
        if(test_2$nearby_beg_new ==0 | test_2$nearby_end_new == 0){
          break
        }
        
        test_3 <- test3(main_brp = main_brp, df_data = df_data,
                       test_1 = test_1, dist_noise = test_2$dist_noise,
                       main_beg = test_2$main_beg_new, 
                       main_end = test_2$main_end_new,
                       nearby_beg = test_2$nearby_beg_new,
                       nearby_end = test_2$nearby_end_new)
        
        out_ind = data.frame(main = main_st, brp = main_brp, nearby = nearby_st,
                             coincident = test_1)
        out_ind_all = out_ind %>% 
          bind_cols(as.data.frame(test_2)) %>% 
          bind_cols(as.data.frame(test_3))
        out <- out %>% 
          rbind(out_ind_all)
      }
    }
  }
}
# ERROR i = 54


# CHARACTERIZATION --------------------------------------------------------


