#' this prog is only used for verifying results with Olivier 

column_classes <- c("character", "Date", "character", rep("numeric",3),
                    rep("Date", 4), rep("numeric", 12))
infor_all = read.table(file = paste0(path_results, "pre_info_test.txt"), 
                       header = TRUE, colClasses = column_classes)
rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE) 
rpt_data$t_break = as.Date(rpt_data$t_break)

nbcsv_min = 249
infor_sel <- infor_all %>% 
  filter(noise < 2, 
         n_main_bef > nbcsv_min, 
         n_nearby_bef > nbcsv_min, 
         n_joint_bef > nbcsv_min, 
         n_main_aft > nbcsv_min, 
         n_nearby_aft > nbcsv_min, 
         n_joint_aft > nbcsv_min)

test <- infor_sel %>% 
  left_join(rpt_data, 
            by = join_by(main == name_main, 
                         brp == t_break,
                         nearby == name_nearby))