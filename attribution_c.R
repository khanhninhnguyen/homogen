
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
#' modify the length of series when it is the begining of cluster:
infor_all$n_main_aft <- infor_all$n_main_aft - infor_all$dist_noise
infor_all$n_joint_aft <- infor_all$n_joint_aft - infor_all$dist_noise

#' add distance & remove the cluster by the first point
rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE) 
distance_list = unique(rpt_data[,c(1,3,13,14)])
infor_all <- infor_all %>% 
  left_join(distance_list, 
            by = join_by(main == name_main, 
                         nearby == name_nearby)) 

#' Print general information 
noise = aggregate(noise ~ main, data = unique(infor_all[,c(1,2,6)]), 
                  FUN = function(x) length(which(x>1)))
total_detection <- aggregate(begin ~ name, data = date_mean, FUN = length) %>%
  filter(name %in% infor_all$main == TRUE) %>%
  mutate(tot = begin - 1) %>%
  left_join(noise, by = join_by(name == main)) %>%
  mutate(final_tot = tot - noise)

print(paste0("total number of detection: ", sum(total_detection$tot)))
print(paste0("total noise are removed: ", sum(total_detection$noise)))
print(paste0("total number of detection after noise removal: ", sum(total_detection$final_tot)))


infor_all <- infor_all %>% filter(noise<2)


# CHARACTERIZATION --------------------------------------------------------

df_infor = list_longest_segment(path_data = path_data_NGL, date_mean = date_mean,
                                path_results = path_results)

characterize(list_infor = df_infor, path_data = path_data_NGL,
             path_results = path_results)


# FGLS TEST ---------------------------------------------------------------
#select first nearby to test by length and distance ---------------------
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
# infor_all <- infor_all %>% 
#   left_join(list_longest_segments,
#             by = join_by(main == main, 
#                         nearby == nearby)) 

#' investigate the length and distance for each changepoint
#' 
#' 
n_min = 200

agg_infor <- infor_all %>%
  mutate(n_joint_min = pmin(n_joint_aft, n_joint_bef)) %>%
  filter(n_main_bef > n_min & n_main_aft > n_min & n_joint_min > n_min) %>% 
  # filter(n_main_bef > n_min & n_main_aft > n_min & dd<50) %>%
  group_by(main, brp) %>%
  summarise(
    mean_n_joint_bef = mean(n_joint_bef, na.rm = TRUE),
    mean_n_joint_aft = mean(n_joint_aft, na.rm = TRUE), 
    mean_d = mean(dd, na.rm = TRUE), 
    total_rows = n()
  )

tested <- aggregate(cbind(UniqueDates = brp) ~ main,
                    data = agg_infor, FUN = function(x) length(unique(x)))

rate = left_join(total_detection, tested, by = join_by(name == main)) %>%
  mutate(r = UniqueDates/final_tot)

print(paste0("total number of fully homogenized: ", length(which(rate$r ==1))))
summary(agg_infor$total_rows)

hist(agg_infor$mean_d, breaks=seq(from=0, to=200, by=10),
     main = "Mean of distance to nearby", 
     xlab = "", 
     ylab = "Frequency")

#' last decision: n_main >200, 10 nearby: d<100 if more than 10, then longest 


# Test for selected  ---------------------
n_min = 200
selected_cases <- infor_all %>%
  mutate(n_joint_min = pmin(n_joint_aft, n_joint_bef)) %>%
  filter(n_main_bef > n_min & n_main_aft > n_min & n_joint_min > n_min) %>% 
  group_by(main, brp) %>%
  group_modify(~select_rows_based_on_conditions(.x)) %>%
  ungroup()
  
write.table(selected_cases, 
            file = paste0(path_results,"list_selected_nmin200_10nearby.txt"),
            row.names = FALSE, quote = FALSE)
#' noise model from data characterization
noise_model_all = read.table(file = paste0(path_results, "order_arma.txt"), 
                             header = FALSE, skip = 1)
noise_model_all[,c(1:3)] <- noise_model_all[,c(1:3)] %>% 
  fill(everything(), .direction = "down")
#" list of segment corresponding to the noise models 
column_classes <- c(rep("character",2), rep("numeric",3),
                    rep("Date", 6))
list_characteried_segments = read.table(file = paste0(path_results, 
                                                  "List_longest_segment.txt"), 
                                    header = TRUE, colClasses = column_classes)
for (i in c(1:2000)) {
  fit.i = list()
  
  main_st = selected_cases$main[i]
  brp = selected_cases$brp[i]
  nearby_st = selected_cases$nearby[i]
  noise.flagged = selected_cases$noise[i]
  dist_noise = selected_cases$dist_noise[i]
  df_data = read_data_new(path_data = path_data_NGL,
                          main_st = main_st, 
                          nearby_st = nearby_st,
                          name_six_diff = name_six_diff)
  
  if(main_st == selected_cases$main[i+1] & brp == selected_cases$brp[i+1]){
    list_ind = c(2:6)
  }else{
    list_ind = c(1:6)
  }
  
  six_noise_models = unlist(noise_model_all[which(
    list_characteried_segments$main == main_st &
      list_characteried_segments$nearby == nearby_st,
  ),])
  
  for (j in list_ind) {
    if(j == 1){
      beg = selected_cases$main_beg_new[i]
      end = selected_cases$main_end_new[i]
    }else if (j == 5){
      beg = selected_cases$nearby_beg_new[i]
      end = selected_cases$nearby_end_new[i]
    }else{
      beg = max(selected_cases$nearby_beg_new[i], selected_cases$nearby_beg_new[i])
      end = min(selected_cases$nearby_end_new[i], selected_cases$nearby_end_new[i])
    }
    
    name_series = name_six_diff[j]
    df = df_data[,c("Date", name_series)] %>% 
      filter(Date >= beg & Date <= end) 
    
    # Split the data frame
    df_before <- df[df$Date <= brp, ]
    df_after <- df[df$Date > brp, ]
    
    # Count non-missing values and limit to first 1000 if necessary
    if (sum(!is.na(df_before[[name_series]])) > 1000) {
      df_before <- df_before %>% filter(Date > (brp - 1000))
    }
    
    if (sum(!is.na(df_after[[name_series]])) > 1000) {
      if (noise.flagged !=0 & j!=5){
        df_after <- df_after %>%
          filter(Date <= (brp + dist_noise + 1 + 1000)) %>%
          filter(Date >= (brp + dist_noise))
      }else{
        df_after <- df_after %>% filter(Date <= (brp + 1000))
      }
    }
    
    new_df <- rbind(df_before, df_after)
    
    df_test <- tidyr::complete(new_df,
                               Date = seq(min(new_df$Date),
                                          max(new_df$Date), 
                                          by = "day"))
    
    
    ind_brp = which(df_test[["Date"]] == brp)
    Data_mod = construct_design(data_df = df_test, 
                                name_series = name_series,
                                break_ind = ind_brp, 
                                one_year = 365)
    
    noise_model = six_noise_models[(j*3-2):(j*3)]
    fit_fgls = FGLS1(design.m = Data_mod,
                     tol= 0.01, 
                     day.list = df_test$Date, 
                     noise.model = noise_model,
                     length.wind0 = 60)
    fit.i[[name_series]] = fit_fgls
  }
  print(i)
  save(fit.i, file = paste0(path_results, main_st, brp, nearby_st, "fgls.RData"))
}



