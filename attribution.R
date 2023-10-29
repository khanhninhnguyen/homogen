#' 5 functions in change-points attribution.
#' INPUT:
#' * 6 series of difference, 
#' * segmentation results of main and nearby station
#'  (date_mean file from segmentation.R) 
#' 
path_data_NGL = "/home/knguyen/data/Data/NGL-attribution/NGL/"
list_six_diff = c("GPS-ERA", "GPS-GPS'", "GPS-ERA'", "ERA-ERA'", "GPS'-ERA'", "GPS'-ERA")
name_six_diff = c("GPS_ERA", "GPS_GPS1", "GPS_ERA1", "ERA_ERA1", "GPS1_ERA1", "GPS1_ERA")
source(file = paste0(path_code, "support_fc.R"))
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
#' 
date_mean_name = paste0("/home/knguyen/Documents/PhD/Results/", "NGL6048.BM_BJ.date_mean.txt")
date_mean = read.table(file = date_mean_name, header = TRUE)
list_brp = extract_list_brp(date_mean) 
list_brp <- list_brp %>%
  group_by(name) %>%
  mutate(cluster_index = get_clusters(brp, 80)) 
infor_all = extract_info_nearby(path_data = path_data,
                                list_brp = list_brp,
                                path_results = path_results)
infor_all = read.table(file = paste0(path_results, "pre_info_test.txt"), 
                       header = TRUE)
# add distance infor 
rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE)
joint_data <- infor_all %>% 
  left_join(rpt_data, 
            by = join_by(main == name_main, 
                         brp == t_break,
                         nearby == name_nearby))
#ERROR
a = full_join(infor_all, rpt_data, by = join_by(main == name_main, 
                                                brp == t_break,
                                                nearby == name_nearby))

check_selected <- function(list_brp, infor_all, nb_min, distance, rate_consecutive){
  filtered <- infor_all %>%
    # Remove noise 
    filter(noise < 2) %>% 
    filter()
  
  
}

df_long <- selected_brp[,c(11:22)] %>%
  gather(key = "variable", value = "value")

# Plot histograms for each column
ggplot(df_long, aes(x = value)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Histograms of Multiple Columns", x = "Value", y = "Frequency")

# CHARACTERIZATION --------------------------------------------------------


