#' this prog is only used for verifying results with Olivier 

column_classes <- c("character", "Date", "character", rep("numeric",3),
                    rep("Date", 4), rep("numeric", 12))
infor_all = read.table(file = paste0(path_results, "pre_info_test.txt"), 
                       header = TRUE, colClasses = column_classes)
rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE) 