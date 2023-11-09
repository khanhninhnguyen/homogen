#' 5 functions in change-points attribution.
#' INPUT:
#' * 6 series of difference, 
#' * segmentation results of main and nearby station
#'  (date_mean file from segmentation.R) 
#' 
path_data_NGL = "/home/knguyen/data/Data/NGL-attribution/NGL/"
list_six_diff = c("GPS-ERA", "GPS-GPS'", "GPS-ERA'", "ERA-ERA'", "GPS'-ERA'", "GPS'-ERA")
reoder_list_name = c("G-E", "G'-E'", "G-G'", "E-E'","G-E'","G'-E")
list_name_test = c("G-E", "G-G'", "G-E'", "E-E'", "G'-E'","G'-E")
name_six_diff = c("GPS_ERA", "GPS_GPS1", "GPS_ERA1", "ERA_ERA1", "GPS1_ERA1", "GPS1_ERA")
source(file = paste0(path_code, "support_fc.R"))
source(file = paste0(path_code, "support_data_characterization.R"))

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
# infor_all = extract_info_nearby(path_data = path_data_NGL,
#                                 list_brp = list_brp,
#                                 path_results = path_results)
column_classes <- c("character", "Date", "character", rep("numeric",3),
                    rep("Date", 4), rep("numeric", 12))
infor_all = read.table(file = paste0(path_results, "pre_info_test.txt"), 
                       header = TRUE, colClasses = column_classes)
# add distance infor, note that the number of points here is computed ----
# by distance between 2 dates, not remove the gaps,
# brps from Olivier is the end of segment, mine is the begining of segment 
rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE) 
distance_list = unique(rpt_data[,c(1,3,13,14)])

infor_all <- infor_all %>% 
  left_join(distance_list, 
            by = join_by(main == name_main, 
                         nearby == name_nearby))

# check selection conditions
check_selected(list_brp = list_brp, infor_all = infor_all,
               nbcsv_min = 200, distance = 200, rate_consecutive = 0.)

# Plot histograms for each column
nbcsv_min = -1
infor_sel <- infor_all %>% 
  filter(noise < 2, 
         n_main_bef > nbcsv_min, 
         n_nearby_bef > nbcsv_min, 
         n_joint_bef > nbcsv_min, 
         n_main_aft > nbcsv_min, 
         n_nearby_aft > nbcsv_min, 
         n_joint_aft > nbcsv_min)
df_long <- infor_all[,c(17:22)] %>%
  gather(key = "variable", value = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Histograms of Multiple Columns", x = "Value", y = "Frequency")

# CHARACTERIZATION -------------------------------------------------------
# only longest segments for each main, nearby and joint series 
df_infor = list_longest_segment(path_data = path_data_NGL, date_mean = date_mean,
                         path_results = path_results)

characterize(list_infor = df_infor, path_data = path_data_NGL,
             path_results = path_results)

# includes 2 steps: 1.fit IGLS to get WLS residual 
# # 2. fit normalized residual to ARMA models 
# nbcsv_min = 200 
# infor_selected <- infor_all %>%
#   filter(noise < 2, 
#          n_main_bef > nbcsv_min, 
#          n_nearby_bef > nbcsv_min, 
#          n_joint_bef > nbcsv_min, 
#          n_main_aft > nbcsv_min, 
#          n_nearby_aft > nbcsv_min, 
#          n_joint_aft > nbcsv_min)
# # keep only 10 closest nearby stations 
# infor_selected <- infor_selected %>%
#   group_by(main, brp) %>%
#   arrange(dd, .by_group = TRUE) %>%
#   slice_head(n = 10) %>%
#   ungroup()

# extract result of data characterization ---------------------------------
column_classes <- c(rep("character",2), rep("numeric",3),
                    rep("Date", 6))
list_selected_segments = read.table(file = paste0(path_results, 
                                                  "List_longest_segment.txt"), 
                       header = TRUE, colClasses = column_classes)
list_selected_segments = list_selected_segments %>% 
  mutate(min_length = apply(list_selected_segments[,c(3:5)], 1, min)) %>% 
  filter(min_length>1000)
# variance 
read_var <- function(path, name_main, name_nearby){
  name_file = paste0("Res_IWLS_", name_main,".", name_nearby, ".RData")
  data_ind = get(load(paste0(path, name_file)))
}
  
  







