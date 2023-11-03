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
infor_all = read.table(file = paste0(path_results, "clean_version/pre_info_test.txt"), 
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

# CHARACTERIZATION --------------------------------------------------------
# includes 2 steps: 1.fit IGLS to get WLS residual 
# 2. fit normalized residual to ARMA models 
nbcsv_min = 200 
infor_selected <- infor_all %>%
  filter(noise < 2, 
         n_main_bef > nbcsv_min, 
         n_nearby_bef > nbcsv_min, 
         n_joint_bef > nbcsv_min, 
         n_main_aft > nbcsv_min, 
         n_nearby_aft > nbcsv_min, 
         n_joint_aft > nbcsv_min)
# keep only 10 closest nearby stations 
infor_selected <- infor_selected %>%
  group_by(main, brp) %>%
  arrange(dd, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# CHECK THIS, NOT TRUE
a = list_longest_segment(path_data = path_data_NGL, date_mean = date_mean,
                         path_results = path_results)
characterize <- function(list_infor, path_data, path_results){
  infor_selected = list_infor
  Four_coef = data.frame(matrix(NA, ncol = 60, nrow = nrow(infor_selected)))
  ARMA_order <- data.frame(matrix(NA, ncol = 18, nrow = nrow(infor_selected)))
  ARMA_coef <- data.frame(matrix(NA, ncol = 24, nrow = nrow(infor_selected)))
  Res_iwls <- list()
  
  for (i in c(1:nrow(infor_selected))) {
    main_st = infor_selected$main[i]
    nearby_st = infor_selected$nearby[i]
    df_data = read_data_new(path_data = path_data_NGL,
                            main_st = main_st, 
                            nearby_st = nearby_st,
                            name_six_diff = name_six_diff)
    df_data$Date <- as.Date(df_data$Date, format = "%Y-%m-%d")
    
    beg <- min(infor_selected$main_beg_new[i], infor_selected$nearby_beg_new[i])
    end <- max(infor_selected$main_end_new[i], infor_selected$nearby_end_new[i])
    
    beg_m <- max(infor_selected$main_beg_new[i], infor_selected$nearby_beg_new[i])
    end_m <- min(infor_selected$main_end_new[i], infor_selected$nearby_end_new[i])
    # replace outside values by NA
    df_data <- df_data %>%
      filter(Date >= beg & Date <= end) %>%
      mutate(
        GPS_ERA = ifelse(Date < infor_selected$main_beg_new[i] |
                           Date > infor_selected$main_end_new[i], NA, GPS_ERA),
        GPS_GPS1 = ifelse(Date < beg_m | Date > end_m, NA, GPS_GPS1),
        GPS_ERA1 = ifelse(Date < beg_m | Date > end_m, NA, GPS_ERA1),
        ERA_ERA1 = ifelse(Date < beg_m | Date > end_m, NA, ERA_ERA1),
        GPS1_ERA1 = ifelse(Date < infor_selected$nearby_beg_new[i] |
                             Date > infor_selected$nearby_end_new[i], NA, GPS1_ERA1),
        GPS1_ERA = ifelse(Date < beg_m | Date > end_m, NA, GPS1_ERA),
      )
    
    Res_IWLS = df_data %>% select(Date)
    brp_ind = which(df_data$Date == infor_selected$brp[i])
    
    for (j in c(1:6)) {
      name.series0 = name_six_diff[j]
      m = construct_design(df_data, name.series = name.series0, break.ind = brp_ind)
      nna_ind = which(!is.na(m$signal))
      
      norm_res <- rep(NA, nrow(m))
      fit_igls_var <- rep(NA, nrow(m))
      fit_igls_res <- rep(NA, nrow(m))
      
      tol0 = 0.01
      if(i == 49 & j ==5){ tol0 = 0.0001 }
      fit_igls = IGLS(design.m = m[nna_ind,], tol =  tol0, day.list = df_data$Date[nna_ind])
      norm_res[nna_ind] = unlist(fit_igls$residual)/sqrt(unlist(fit_igls$var))
      arima_fit = fit.arima(norm_res)
      
      fit_igls_var[nna_ind] <- unlist(fit_igls$var)
      fit_igls_res[nna_ind] <- unlist(fit_igls$residual)
      
      Res_IWLS[, paste0(name.series0, '_var')] <- fit_igls_var
      Res_IWLS[, paste0(name.series0, '_res')] <- fit_igls_res
      
      ARMA_order[i,c((3*j-2):(3*j))] = arima_fit$pq
      ARMA_coef[i,c((4*j-3):(4*j))] = round(arima_fit$coef, digits = 4)
      Four_coef[i,c((10*j-9):(10*j))] = round(fit_igls$coefficients, digits = 4)
    }
    # Res_iwls[[paste(
    #   infor_all[i, 1], 
    #   format(infor_all[i, 2], "%Y-%m-%d"), 
    #   infor_all[i, 3], sep = ".")]] <- Res_IWLS
    name_case = paste(infor_all[i, 1],
                      format(infor_all[i, 2], "%Y-%m-%d"), 
                      infor_all[i, 3], sep = ".")
    save(Res_IWLS, 
         file = paste0(path_results, "Res_IWLS_", name_case, ".RData"))
    
    print(i)
  }
  
  write.table(ARMA_order, file = paste0(path_results, "order_arma.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(ARMA_coef, file = paste0(path_results, "coef_arma.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(Four_coef, file = paste0(path_results, "Four_coef.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  #' Modify a bit the final result
  order_arma = read.table(file = paste0(path_results, "order_arma.txt"),
                          header = TRUE)
  coef_arma = read.table(file = paste0(path_results, "coef_arma.txt"),
                         header = TRUE)
  
  order_arma_m = order_arma[,-seq(2,18,3)]
  colnames(order_arma_m) = c(rbind(outer(c("AR-", "MA-"), 
                                         name_six_diff, paste0)))
  coef_arma_m = coef_arma[,-seq(2,24,2)]
  colnames(coef_arma_m) = c(rbind(outer(c("Phi-", "Theta-"), 
                                        name_six_diff, paste0)))
  
  write.table(order_arma_m, file = paste0(path_results, "Order_ARMA.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(coef_arma_m, file = paste0(path_results, "Coeff_ARMA.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# extract result of data characterization ---------------------------------
# NOTE PLOT AND NUMBER SHOULD REPORT ONLY 1 FOR EACH PAIR OF MAIN-NEARBY 
# --> CAN CHECK THE CONSISTENCY IN NOISE MODEL IN EACH PAIR 

  
  




