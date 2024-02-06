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

SD_infor = get(load(file = paste0(path_results, "mean_range_SD.RData")))
list_SD_infor = list_selected_segments %>%
  filter(length_joint>1000)

mean_sd = SD_infor$mean_sd
range_sd = SD_infor$range_sd

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


# old_list = get(load(file = paste0(path_results, "list.segments.selected.RData")))
# ind_coin = which(list_selected_segments$nearby %in% old_list$main)
# var_list = mean_sd[ind_coin,]
# df = cbind(list_selected_segments[ind_coin,], var_list)
# df$duplicated_marker <- duplicated(df$main)
# df$main[df$duplicated_marker] <- NA

# plot sd - distance 

df_all = cbind(list_SD_infor, mean_sd)

main_duplicated = which(duplicated(df_all$main))
nearby_duplicated = which(duplicated(df_all$nearby))

df_all$GPS_ERA[main_duplicated] <- NA
df_all$GPS1_ERA1[nearby_duplicated] <- NA

df_plot = df_all[,c("dd", name_six_diff)]

df_long <- df_plot %>%
  pivot_longer(cols = 2:7, names_to = "Variable", values_to = "Value")

# Plot with ggplot2

# Create a named vector for subtitles
subtitles <- c("GPS_ERA" = paste0("GPS-ERA (",
                                  length(df_plot$GPS_ERA[!is.na(df_plot$GPS_ERA)]),
                                  "cases)"), 
               "GPS_GPS1" = paste0("GPS-GPS' (",
                                   length(df_plot$GPS_GPS1[!is.na(df_plot$GPS_GPS1)]),
                                   "cases)"), 
               "GPS_ERA1" = paste0("GPS-ERA' (", 
                                   length(df_plot$GPS_ERA1[!is.na(df_plot$GPS_ERA1)]),
                                   "cases)"), 
               "ERA_ERA1" = paste0("ERA-ERA' (", 
                                   length(df_plot$ERA_ERA1[!is.na(df_plot$ERA_ERA1)]),
                                   "cases)"), 
               "GPS1_ERA1" = paste0("GPS'-ERA' (", 
                                    length(df_plot$GPS1_ERA1[!is.na(df_plot$GPS1_ERA1)]),
                                    "cases)"), 
               "GPS1_ERA" = paste0("GPS'-ERA  (", 
                                   length(df_plot$GPS1_ERA[!is.na(df_plot$GPS1_ERA)]),
                                   "cases)"))

p <- ggplot(df_long, aes(x = dd, y = Value)) +
  geom_point(size = 0.3) +
  facet_wrap(~ Variable, scales = "free_x", labeller = labeller(Variable = subtitles))+
  theme_minimal() +
  labs(x = "Distance", y = "Value")
ggsave(paste0(path_results,"attribution/SD_distance.jpg" ), 
       plot = p, width = 22, height = 16, units = "cm", dpi = 1200)

# plot length distribution 

p <- df_all[,c(3:5)] %>%
  reshape2::melt() %>%
  ggplot( aes(x=value)) + theme_bw()+
  geom_histogram( position = 'identity', bins = 100) +
  facet_wrap(~ variable, scales = "free_x")
ggsave(paste0(path_results,"attribution/length_histogram.jpg" ), 
       plot = p, width = 14, height = 5, units = "cm", dpi = 1200)
# arma models -------------------------------------------------

text1 = "Distance < 50 km"
text2 = "Distance >= 50 km"

order_arma = read.table(file = paste0(path_results, "order_arma.txt"),
                        header = TRUE,
                        col.names = as.vector(outer(c("p","d","q"),
                                                    name_six_diff,
                                                    paste, sep=".")))

# make sure that filtered dataframe has enough data for G-E 
order_arma[,c(1:3)] <- order_arma[,c(1:3)] %>%
  fill(everything(), .direction = "down")

ind_selected = which(list_selected_segments$dd < 50)

infor_selected = list_selected_segments[ind_selected,]
order_arma_filtered = order_arma[ind_selected,]
order_arma_filtered[which(duplicated(infor_selected$main)),
                    paste(c("p","d","q"), name_six_diff[1], sep = ".")] <- NA
order_arma_filtered[which(duplicated(infor_selected$nearby)),
                    paste(c("p","d","q"),name_six_diff[5], sep = ".")] <- NA

plot_dist_model(arma_order = order_arma_filtered, 
                name_fig = text1, 
                name_six_diff, sub_title = text1, tag_fig = "(a)" )

# arma coefficients -------------------------------------------
coef_arma = read.table(file = paste0(path_results, "coef_arma.txt"), 
                       header = TRUE, 
                       col.names = as.vector(outer(c("Phi1", "Phi2","Theta1","Theta2"),
                                                   name_six_diff,
                                                   paste, sep=".")))

# make sure that filtered dataframe has enough data for G-E 
coef_arma[,c(1:4)] <- coef_arma[,c(1:4)] %>%
  fill(everything(), .direction = "down")

ind_selected = which(list_selected_segments$dd >= 50)

infor_selected = list_selected_segments[ind_selected,]
order_arma_filtered = order_arma[ind_selected,]
order_arma_filtered[which(duplicated(infor_selected$main)), c(1:3)] <- NA
order_arma_filtered[which(duplicated(infor_selected$nearby)), c(13:15)] <- NA

coef_arma_filtered = coef_arma[ind_selected,]
coef_arma_filtered[which(duplicated(infor_selected$main)), 
                   paste(c("Phi1", "Phi2","Theta1","Theta2"), name_six_diff[1], sep = ".")] <- NA
coef_arma_filtered[which(duplicated(infor_selected$nearby)), 
                   paste(c("Phi1", "Phi2","Theta1","Theta2"), name_six_diff[5], sep = ".")] <- NA



plot_arma_coef(arma_order = order_arma_filtered, 
               arma_coef = coef_arma_filtered, 
               list_infor = infor_selected,
               name_fig = text2,
               sub_title = text2,
               tag_fig = "(b)")
# Plot theta as a function of phi from ARMA(1,1)
arma_order = order_arma[which(list_selected_segments$length_joint>1000),]
arma_coef = coef_arma[which(list_selected_segments$length_joint>1000),seq(1,24,2)]

list_model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)")

length_data = nrow(arma_order)
six_model = data.frame(matrix(NA, ncol = 6, nrow = length_data))
for (i in 1:6) {
  six_model[,i] = sapply(c(1:length_data), function(x) {
    model.iden(as.numeric(unlist(arma_order[x, (3*i-2):(3*i)])))
  })
}
colnames(six_model) <- name_six_diff
a = list_selected_segments[which(list_selected_segments$length_joint>1000),]
six_model[which(duplicated(a$main)),1] <- NA
six_model[which(duplicated(a$nearby)),1] <- NA

for (i in c(1:6)) {
  df = arma_coef[,(2*i-1):(2*i)]
  df <- df[which(six_model[,i] == "ARMA(1,1)"),]
  p <- ggplot(data = df, aes_string(x=names(df)[1], y=names(df)[2]))+
    theme_bw()+
    geom_point(size = 0.2, aes(col = as.factor(class)))+
    xlim(-1, 1) + ylim(-1, 1) +
    geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "blue")
  ggsave(paste0(path_results,"attribution/coeff_ARMA",name_six_diff[i],".jpg" ), 
         plot = p, width = 10, height = 8, units = "cm", dpi = 600)
}

df_sd = cbind(list_SD_infor, mean_sd)
df_all = left_join(list_selected_segments, df_sd, by = c("main", "nearby"))

df$class = 0
df$class[which(df$Theta1.GPS1_ERA>0 & df$Theta1.GPS1_ERA<0.35)] <- 1
# df$class[which(abs(df$Theta1.GPS1_ERA)>0.7 | abs(df$Theta1.GPS1_ERA>0.7))] <- 1
df$l = list_selected_segments$dh[which(six_model[,i] == "ARMA(1,1)")]
df$sd = df_all$GPS1_ERA[which(six_model[,i] == "ARMA(1,1)")]


df %>%
  ggplot( aes(x=l, fill=as.factor(class))) + theme_bw()+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
  scale_fill_manual(values=c("#69b3a2", "#404080"))

df %>%
  ggplot(aes(x = l, y = sd, col = as.factor(class))) +
  theme_bw() +
  geom_point()

#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' test results
#' 



# test results ------------------------------------------------------------
library(gridExtra)
path_restest <- paste0(path_results,"attribution/predictive_rule/")
List_main = read.table(file = paste0(path_restest,"list_selected_nmin200_10nearby.txt"), 
                       header = TRUE, 
                       stringsAsFactors = FALSE) 
name_six_diff = c("GPS_ERA", "GPS_GPS1", "GPS_ERA1", "ERA_ERA1", "GPS1_ERA1", "GPS1_ERA")

name.version ="FGLS_jump_tvalue.txt"
name.results <- paste0(path_restest, name.version) # name of test result file
Data_Res_Test0 <- read.table(name.results,
                             header = TRUE,
                             stringsAsFactors = FALSE) 

for(suffix in name_six_diff) {
  jump_col <- paste0("Jump_", suffix)
  tvalue_col <- paste0("Tvalue_", suffix)
  
  # Create new ratio column
  Data_Res_Test0[[paste0("Std.err_", suffix)]] <- Data_Res_Test0[[jump_col]] / Data_Res_Test0[[tvalue_col]]
}

Data_Res_Test <- cbind(List_main[,c("main", "brp", "nearby")],
                       Data_Res_Test0[,1:6]) %>%
  mutate(brp = as.Date(List_main$brp, format="%Y-%m-%d")) 



suspect = Data_Res_Test0 %>%
  filter(abs(Jump_GPS_GPS1)<0.01)
suspect_info = List_main[which(abs(Data_Res_Test0$Jump_GPS_GPS1)<0.01),]
# distribution -------------------

# Transform the dataset to a long format for easier plotting with ggplot2
# df_plot = Data_Res_Test[which(List_main$dd>50),]
Data_Res_Test_fillNA <- Data_Res_Test %>%
  arrange(main, brp) %>%
  group_by(main, brp) %>%
  mutate(Std.err_GPS_ERA = if_else(is.na(Jump_GPS_ERA), 
                                  lag(Jump_GPS_ERA, order_by = brp, default = NA), 
                                  Jump_GPS_ERA)) %>%
  fill(Jump_GPS_ERA, .direction = "downup") %>%
  ungroup()

df_plot = Data_Res_Test %>% 
  filter(Jump_GPS_ERA<0)

# df_long <- pivot_longer(df_plot, cols = starts_with("Tvalue"), names_to = "Variable", values_to = "Value")
df_long <- pivot_longer(df_plot, cols = starts_with("Jump"), names_to = "Variable", values_to = "Value")

# Create the plots
ggplot(df_long, aes(x = Value)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") + # Adjust binwidth as needed
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  labs(title = "Distribution of Std.err ( d>25 km)", x = "Value", y = "Frequency") +
  # labs(title = "Distribution of T-value Columns", x = "Value", y = "Density") +
  theme_minimal() +
  xlim(0,0.7)
# Inspect the suspected cases ----------------------------
# 
# suspect1 = which(abs(Data_Res_Test_fillNA$Jump_GPS_ERA)>3)
# suspect2 = which(abs(Data_Res_Test_fillNA$Jump_GPS_GPS1)>3)
# suspect3 = which(abs(Data_Res_Test_fillNA$Jump_GPS_ERA1)>3)

suspect_case = unique(suspect1, suspect2, suspect3)
suspect_case =  which(apply(Data_Res_Test[,4:9], 1, function(x) any(x > 30)))

for (i in suspect_case) {
  main_st = Data_Res_Test$main[i] 
  brp_test = Data_Res_Test$brp[i]
  name_nearby_full = Data_Res_Test %>% 
    filter(main == main_st,
           brp == brp_test) %>%
    select(nearby)  %>%
    pull %>%
    tail(n = 1)
  
  plot_test_res(main_st = Data_Res_Test$main[i] ,
                brp = Data_Res_Test$brp[i], 
                nearby_st = Data_Res_Test$nearby[i],  
                main_beg = List_main$main_beg_new[i],
                main_end = List_main$main_end_new[i], 
                nearby_beg = List_main$nearby_beg_new[i], 
                nearby_end = List_main$nearby_end_new[i],
                name_nearby_full = name_nearby_full,
                name_six_diff,
                path_data_NGL)
  plot_full_series(main_st = Data_Res_Test$main[i], 
                   nearby_st = Data_Res_Test$nearby[i], 
                   path_data_NGL = path_data_NGL, 
                   date_mean = date_mean,
                   name_six_diff = name_six_diff)
}


