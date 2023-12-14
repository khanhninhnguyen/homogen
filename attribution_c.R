
# clear version -----------------------------------------------------------

#' 5 functions in change-points attribution.
#' INPUT:
#' * 6 series of difference, 
#' * segmentation results of main and nearby station
#'  (date_mean file from segmentation.R) 
#' procedure:
#' 1. extract list of cases can be tested (all information)
#' 2. characterization on the longest segment (min = 1000) in each pair 
#' main-nearby 
#' 3. using distance and noise (estimated from the 2nd step) to select the case 
#' used for the FGLS test 
#' 
#' 
# Read data and functions ------------------------------------------------
path_data_NGL = "/home/knguyen/data/Data/NGL-attribution/NGL/"
list_six_diff = c("GPS-ERA", "GPS-GPS'", "GPS-ERA'", "ERA-ERA'", "GPS'-ERA'", "GPS'-ERA")
reoder_list_name = c("G-E", "G'-E'", "G-G'", "E-E'","G-E'","G'-E")
list_name_test = c("G-E", "G-G'", "G-E'", "E-E'", "G'-E'","G'-E")
name_six_diff = c("GPS_ERA", "GPS_GPS1", "GPS_ERA1", "ERA_ERA1", "GPS1_ERA1", "GPS1_ERA")
source(file = paste0(path_code, "support_fc.R"))
# source(file = paste0(path_code, "support_data_characterization.R"))
source(file = paste0(path_code, "support_FGLS.R"))
# order: GPS-ERA, GPS-GPS', GPS-ERA', ERA-ERA', GPS'-ERA', GPS'-ERA

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

#' uncomment this when run the first time, result stored in pre_infor_test.txt
#' infor_all = extract_info_nearby(path_data = path_data_NGL,
#'                                 list_brp = list_brp,
#'                                 path_results = path_results)
column_classes <- c("character", "Date", "character", rep("numeric",3),
                    rep("Date", 4), rep("numeric", 12))
infor_all = read.table(file = paste0(path_results, "pre_info_test.txt"), 
                       header = TRUE, colClasses = column_classes)

#' add distance 
rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE) 
distance_list = unique(rpt_data[,c(1,3,13,14)])
infor_all <- infor_all %>% 
  left_join(distance_list, 
            by = join_by(main == name_main, 
                         nearby == name_nearby))







# CHARACTERIZATION --------------------------------------------------------

df_infor = list_longest_segment(path_data = path_data_NGL, date_mean = date_mean,
                                path_results = path_results)

characterize(list_infor = df_infor, path_data = path_data_NGL,
             path_results = path_results)


# FGLS TEST ---------------------------------------------------------------
#' select first nearby to test by SD and distance 
#' read data
SD_infor = get(load(file = paste0(path_results, "mean_range_SD.RData")))

column_classes <- c(rep("character",2), rep("numeric",3),
                    rep("Date", 6))
list_longest_segments = read.table(file = paste0(path_results, 
                                                  "List_longest_segment.txt"), 
                                    header = TRUE, colClasses = column_classes) %>%
  filter(length_joint>1000) %>% 
  bind_cols(set_names(SD_infor$mean_sd, paste0("SD_", names(SD_infor$mean_sd)))) %>%
  select(-c(length_main, length_joint, length_nearby,
            beg_main, beg_nearby, beg_joint, 
            end_main, end_nearby, end_joint))
#' filtering 
#' 
infor_all <- infor_all %>% 
  left_join(list_longest_segments,
            by = join_by(main == main, 
                        nearby == nearby)) 

# if there are more than 10 nearby: select the length and noise 
list_selected <- infor_all %>% 
  filter(n_joint_bef>=100 & n_joint_aft>=100) %>%
  group_by(main, brp) %>%
  filter(if (n() >= 10)) {
    n_joint_bef > 400 & n_joint_aft > 400 & dd < 100
  } else {
    TRUE
  })
  # mutate(n_joint = n_joint_bef + n_joint_aft) %>%
  # group_by(main, brp) %>%
  # filter(if (n() >= 10 && any(n_joint_bef > 400) && any(n_joint_aft > 400)) {
  #   n_joint_bef > 400 & n_joint_aft > 400 & dd < 100
  # } else {
  #   TRUE
  # }) %>%
  # arrange(SD_GPS_GPS1) %>%  # Arrange by 'sd' within each group
  # slice_head(n = min(10, n()))
# to check percentage of fully homogenized stations 
# a = unique(infor_all[,c(1,2,6)])
# noise = aggregate(noise ~ main, data = a, FUN = function(x) length(which(x!=0)))
# all <- aggregate(begin ~ name, data = date_mean, FUN = length)
# tested <- aggregate(cbind(UniqueDates = brp) ~ main,
#                     data = list_selected, FUN = function(x) length(unique(x))) %>%
#   left_join(all,
#             by = join_by(main == name)) %>%
#   left_join(noise,
#             by = join_by(main == main)) %>%
#   mutate(r = UniqueDates/(begin-1-noise))
# table(tested$r == 1)




