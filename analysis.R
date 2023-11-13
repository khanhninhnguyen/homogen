#' ANALYZE RESULT 
#' data characterization
#' 
source(file = paste0(path_code, "support_data_characterization.R"))

column_classes <- c(rep("character",2), rep("numeric",3),
                    rep("Date", 6))
list_selected_segments = read.table(file = paste0(path_results, 
                                                  "List_longest_segment.txt"), 
                                    header = TRUE, colClasses = column_classes)
rpt_data <- read.table(file = paste0(path_data, "support/liste_main20yr_1nearby_200km_500m_np250_nd250.rpt"),
                       header = TRUE, check.names = FALSE) 
distance_list = unique(rpt_data[,c(1,3,13,14)])
list_selected_segments <- list_selected_segments %>% 
  left_join(distance_list, 
            by = join_by(main == name_main, 
                         nearby == name_nearby))

# variance ---------------------------------

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


# arma models -------------------------------------------------

text1 = "Distance < 25 km"
text2 = "Distance > 25 km"

order_arma = read.table(file = paste0(path_results, "order_arma.txt"),
                        header = TRUE)

plot_dist_model(arma_order = order_arma[which(list_selected_segments$dd < 25),], 
                name_fig = text1, 
                name_six_diff, sub_title = text1, tag_fig = "(a)" )

plot_dist_model(arma_order = order_arma[which(list_selected_segments$dd > 25),], 
                name_fig = text2, 
                name_six_diff, sub_title = text2, tag_fig = "(b)" )
# arma coefficients -------------------------------------------
coef_arma = read.table(file = paste0(path_results, "coef_arma.txt"), 
                       header = TRUE)

plot_arma_coef(arma_order = order_arma[which(list_selected_segments$dd < 25),], 
               arma_coef = coef_arma[which(list_selected_segments$dd < 25),], 
               list_infor = list_selected_segments[which(list_selected_segments$dd < 25),],
               name_fig = "test", sub_title = text1, tag_fig = "(a)")
# DEAL WITH NA IN GPS-ERA
#' 
#' 
#' 
#' 
#' 
#' 
#' 