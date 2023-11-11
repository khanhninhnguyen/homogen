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
# source(file = paste0(path_code, "support_data_characterization.R"))
source(file = paste0(path_code, "support_FGLS.R"))
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

rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE) 
distance_list = unique(rpt_data[,c(1,3,13,14)])

# add distance infor, note that the number of points here is computed ----
# by distance between 2 dates, not remove the gaps,
# brps from Olivier is the end of segment, mine is the begining of segment 

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
extract_info_SD <- function(path, list_segment, name_six_diff, path_results){
  n = nrow(list_segment)
  a <- data.frame(matrix(NA, nrow = n, ncol = 6))
  colnames(a) <- name_six_diff
  mean_sd <- a
  range_sd <- a
  
  for (i in c(1:n)) {
    var_6diff = read_var(path = path, 
                         name_main = list_segment$main[i],
                         name_nearby = list_segment$nearby[i],
                         name_six_diff = name_six_diff)
    mean_sd[i,] = colMeans(sqrt(var_6diff[,c(name_six_diff)]),
                            na.rm = TRUE)
    range_sd[i,] = sapply(name_six_diff, function(x){
      range_cal(variable = sqrt(var_6diff[[x]]), 
                day.list = var_6diff$Date)})
    print(i)
  }
  write.table(mean_sd, file = paste0(path_results, "mean_sd.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(range_sd, file = paste0(path_results, "range_sd.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  out = list(mean_sd = mean_sd, range_sd = range_sd)
  return(out)
}  
  
a = extract_info_SD(path = paste0(path_results,"characterization/"),
                list_segment = list_selected_segments,
                name_six_diff = name_six_diff)
mean_var = a$mean_sd
range_var = a$range_sd
distance = list_selected_segments %>% 
  left_join(distance_list, 
            by = join_by(main == name_main, 
                         nearby == name_nearby))
ind1 = which(distance$dd<50)
ind2 = which(distance$dd>=50)

ind = ind1
res = data.frame(mean = colMeans(mean_var[ind,], na.rm = TRUE),
                 sd.mean = apply(mean_var[ind,], 2 , sd, na.rm = TRUE), 
                 range = colMeans(range_var[ind,]/ mean_var[ind,], na.rm = TRUE),
                 sd.range = apply(range_var[ind,]/ mean_var[ind,], 2 , sd, na.rm = TRUE),
                 name = name_six_diff)

ind = ind2
res = data.frame(mean = colMeans(mean_var[ind,], na.rm = TRUE),
                 sd.mean = apply(mean_var[ind,], 2 , sd, na.rm = TRUE), 
                 range = colMeans(range_var[ind,]/ mean_var[ind,], na.rm = TRUE),
                 sd.range = apply(range_var[ind,]/ mean_var[ind,], 2 , sd, na.rm = TRUE),
                 name = name_six_diff)

# write.table(a, file = paste0(path_results,"attribution/sd_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
colMeans(range_var, na.rm = TRUE)
apply(mean_var_la, 2 , sd, na.rm = TRUE)







# select the list of cases to do the test ---------------------------------
# by limiting number of nearby vs distance, noise 
SD_info = get(load(file = paste0(path_results, "mean_range_SD.RData")))
SD_info$mean_sd$main = list_selected_segments$main
SD_info$mean_sd$nearby = list_selected_segments$nearby

nbcsv_min = 250
infor_selected <- infor_all %>%
  filter(noise < 2,
         n_main_bef > nbcsv_min,
         n_nearby_bef > nbcsv_min,
         n_joint_bef > nbcsv_min,
         n_main_aft > nbcsv_min,
         n_nearby_aft > nbcsv_min,
         n_joint_aft > nbcsv_min)

number_nearby <- infor_selected %>%
  group_by(main, brp) %>%
  summarise(count = n())
# investigate ----

infor_selected <- infor_selected %>%
  left_join(distance_list, 
            by = join_by(main == name_main, 
                         nearby == name_nearby)) %>%
  left_join(SD_info$mean_sd,
            by = join_by(main == main,
                          nearby == nearby)) %>%
  left_join(number_nearby,
                by = join_by(main == main,
                             brp == brp)) %>%
  filter(count>10)

data_sel_cri <- infor_selected[,c(1:3,13,16,19,22,23,26,29)]

list_10dd <- data_sel_cri %>%
  group_by(main, brp) %>%
  arrange(dd) %>%
  slice_head(n = 10)  %>%
  mutate(cri = "distance")

list_10noise <- data_sel_cri %>%
  group_by(main, brp) %>%
  arrange(GPS_GPS1) %>%
  slice_head(n = 10) %>%
  mutate(cri = "noise_gg")
list_10noise1 <- data_sel_cri %>%
  group_by(main, brp) %>%
  arrange(GPS1_ERA1) %>%
  slice_head(n = 10) %>%
  mutate(cri = "noise_nb")

data_plot = rbind(list_10dd, list_10noise, list_10noise1)
long_df <- data_plot %>%
  pivot_longer(cols = 4:10, names_to = "variable", values_to = "value")
# visualize
ggplot(long_df, aes(x = variable, y = value, color = cri)) +
  theme_bw() + 
  geom_boxplot() +
  facet_wrap(~variable, scales = "free") +
  labs(fill = "Group")
ggplot(data = data_plot, aes(x = dd, col = as.factor(cri))) + 
  theme_bw() + 
  geom_boxplot(aes(y = GPS_GPS1))

# extract final list ----

final_list <- infor_selected %>%
  left_join(distance_list, 
            by = join_by(main == name_main, 
                         nearby == name_nearby)) %>%
  left_join(SD_info$mean_sd,
            by = join_by(main == main,
                         nearby == nearby)) %>%
  left_join(number_nearby,
            by = join_by(main == main,
                         brp == brp)) %>%
  group_by(main, brp) %>%
  arrange(GPS_GPS1) %>%
  slice_head(n = 10) 

# select europe 

latlon = read.table(file = paste0(path_data, "support/",
                                  "gps_sta_NGL_dates_np_10.0yr_0.10gaps_bkl.txt"),
                    header = TRUE, check.names = FALSE)
latlon$name_ngl <- tolower(latlon$name_ngl)

final_list_location <- final_list %>% 
  left_join(latlon[,c(1:3)], by = join_by(main == name_ngl)) 
final_list_location$lon <- ifelse(final_list_location$lon >0 & final_list_location$lon < 180, 
                                  final_list_location$lon, 
                                  -(360 -final_list_location$lon))  
con = final_list_location$lon > -10 & final_list_location$lon < 50 &
  final_list_location$lat > 30 & final_list_location$lat < 75
final_list_location$select <- ifelse(con, 1, 0)
library(ggplot2)
library(maps)
# Basic world map
world_map <- map_data("world")

# Plot the map and add points
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
  geom_point(data = final_list_location, aes(x = lon, y = lat, color = as.factor(select))) +
  theme_minimal() +
  labs(title = "Map of selected stations")

final_list_location <- final_list_location %>%
  filter(select ==1)
write.table(final_list_location, file = paste0(path_results, "final_list_Europe.txt"), 
            sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# SIGNIFICANCE TEST -------------------------------------------------------
#' We need: infor file, model ARMA file, 
column_classes <- c("character", "Date", "character", rep("numeric",3),
                    rep("Date", 4), rep("numeric", 24))
infor_all = read.table(file = paste0(path_results, "final_list_Europe.txt"), 
                       header = TRUE, colClasses = column_classes)
noise_model_all = read.table(file = paste0(path_results, "order_arma.txt"), 
                             header = FALSE, skip = 1)
noise_model_all[,c(1:3)] <- noise_model_all[,c(1:3)] %>% 
  fill(everything(), .direction = "down")

column_classes <- c(rep("character",2), rep("numeric",3),
                    rep("Date", 6))
list_selected_segments = read.table(file = paste0(path_results, 
                                                  "List_longest_segment.txt"), 
                                    header = TRUE, colClasses = column_classes)
for (i in c(1:1000)) {
  fit.i = list()
  
  main_st = infor_all$main[i]
  brp = infor_all$brp[i]
  nearby_st = infor_all$nearby[i]
  df_data = read_data_new(path_data = path_data_NGL,
                          main_st = main_st, 
                          nearby_st = nearby_st,
                          name_six_diff = name_six_diff)
  
  if(main_st == infor_all$main[i+1] & brp == infor_all$brp[i+1]){
    list_ind = c(2:6)
  }else{
    list_ind = c(1:6)
  }
  
  six_noise_models = unlist(noise_model_all[which(
    list_selected_segments$main == main_st &
      list_selected_segments$nearby == nearby_st,
  ),])
  
  for (j in list_ind) {
    if(j == 1){
      beg = infor_all$main_beg_new[i]
      end = infor_all$main_end_new[i]
    }else if (j == 5){
      beg = infor_all$nearby_beg_new[i]
      end = infor_all$nearby_end_new[i]
    }else{
      beg = max(infor_all$nearby_beg_new[i], infor_all$nearby_beg_new[i])
      end = min(infor_all$nearby_end_new[i], infor_all$nearby_end_new[i])
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
      df_after <- df_after %>% filter(Date <= (brp + 1000))
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
    fit.i[[name_series]] = fit.fgls
  }
  print(i)
  save(fit.i, file = paste0(path_results, main_st, brp, nearby_st, "fgls.RData"))
}

