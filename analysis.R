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
List_main = read.table(file = paste0(path_results,"list_selected_nmin200_10nearby.txt"), 
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

# suspect = Data_Res_Test0 %>%
#   filter(abs(Jump_GPS_GPS1)<0.01)
# suspect_info = List_main[which(abs(Data_Res_Test0$Jump_GPS_GPS1)<0.01),]
# distribution -------------------

# Transform the dataset to a long format for easier plotting with ggplot2
# df_plot = Data_Res_Test[which(List_main$dd>50),]
Data_Res_Test_fillNA <- Data_Res_Test %>%
  arrange(main, brp) %>%
  group_by(main, brp) %>%
  mutate(Std.err_GPS_ERA = if_else(is.na(Tvalue_GPS_ERA), 
                                  lag(Tvalue_GPS_ERA, order_by = brp, default = NA), 
                                  Tvalue_GPS_ERA)) %>%
  fill(Tvalue_GPS_ERA, .direction = "downup") %>%
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

####' Distribution of G-E in 2 sides 

df_plot = data.frame(Raw = na.omit(Data_Res_Test0$Tvalue_GPS_ERA)) %>%
  mutate(Value = abs(Raw))
df_plot$sign <- ifelse(df_plot$Raw < 0, "Negative", "Positive")

df_plot %>%
  
  ggplot( aes(x=Value, fill=sign)) + theme_bw()+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + 
  labs(title = "Histogram of the amplitude of T-values separated by sign")


# population of tests -----------------------------------------------------
df = Data_Res_Test_fillNA[which(fix_case$Fix == 0),-10] 
df <- df %>%
  group_by(main, brp) %>%
  mutate(Tvalue_GPS_ERA = ifelse(duplicated(Tvalue_GPS_ERA), NA, Tvalue_GPS_ERA)) %>%
  ungroup()
# Reshape data to long format
df_long <- df[, -c(1:3)] %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
  drop_na()  # Drop NA values to avoid including them in the categorization

# Categorize each value
df_long$Category <- case_when(
  df_long$Value < -1.96 ~ "-1",
  df_long$Value >= -1.96 & df_long$Value <= 1.96 ~ "0",
  df_long$Value > 1.96 ~ "1",
  TRUE ~ "Other"
)

# Count the occurrences of each category for each variable
df_summary <- df_long %>%
  count(Variable, Category)
df_summary$Variable <- as.character(sub("Tvalue_", "", df_summary$Variable))
df_summary$Variable = factor(df_summary$Variable,  
                            levels = c("GPS_ERA", "GPS_ERA1", "GPS_GPS1", "GPS1_ERA1", "GPS1_ERA", "ERA_ERA1"))
df_summary$Category = factor(df_summary$Category, levels = c("-1", "0", "1"))
# Plot
p <- ggplot(df_summary, aes(fill=Category, y=n, x=Variable)) + 
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  geom_text(aes(label = n),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count") + 
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
ggsave(paste0(path_results,"attribution/pop_significance_level1.jpg" ), 
       plot = p, 
       width = 8,
       height = 5,
       units = "cm",
       dpi = 1200,
       bg="white")


# Inspect the suspected cases ----------------------------
# 
# suspect1 = which(abs(Data_Res_Test_fillNA$Jump_GPS_ERA)>3)
# suspect2 = which(abs(Data_Res_Test_fillNA$Jump_GPS_GPS1)>3)
# suspect3 = which(abs(Data_Res_Test_fillNA$Jump_GPS_ERA1)>3)

suspect_case = unique(suspect1, suspect2, suspect3)
# suspect_case =  which(apply(Data_Res_Test[,4:9], 1, function(x) any(x > 30)))
suspect_case = which(abs(Data_Res_Test0$Tvalue_GPS_ERA) <1.96)
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
                name_nearby_full = name_nearby_full,
                name_six_diff,
                path_FGLS = paste0(path_results, "attribution/FGLS_cor/"),
                path_data_NGL,
                distance = round(List_main$dd[i], digits = 1))
  plot_full_series(main_st = Data_Res_Test$main[i],
                   nearby_st = Data_Res_Test$nearby[i],
                   path_data_NGL = path_data_NGL,
                   date_mean = date_mean,
                   name_six_diff = name_six_diff,
                   distance = round(List_main$dd[i], digits = 1))
}




# comparison result of prediction -----------------------------------------


# Define variable names and read results 
ori_col_name <- "Ori"
ver1_col_name <- "Ver1"
title_name = "Results when R =400 of"

ori = read.table(file = paste0(path_restest, 'original', "/FinalTable.txt"))
ver1 = read.table(file = paste0(path_restest, 'ver1', "/FinalTable.txt"))

## Create a new dataframe with specified column names
df = data.frame(ori = ori$pred.y, ver1 = ver1$pred.y)
names(df) <- c(ori_col_name, ver1_col_name)  # Apply the variable names

all_values <- sort(unique(c(df[[ori_col_name]], df[[ver1_col_name]])))


common_freqs <- df %>%
  filter(!!sym(ori_col_name) == !!sym(ver1_col_name)) %>%
  count(!!sym(ori_col_name)) %>%
  rename(Value = ori_col_name, CommonFreq = n)

## Step 2: Calculate individual frequencies for visualization
df_counts <- df %>%
  pivot_longer(cols = all_of(c(ori_col_name, ver1_col_name)), names_to = "Column", values_to = "Value") %>%
  group_by(Column, Value) %>%
  summarise(Frequency = n(), .groups = 'drop') %>%
  ungroup()

## Join the common frequencies with the df_counts for plotting
df_counts <- df_counts %>%
  left_join(common_freqs, by = "Value")

## Step 3: Plotting
ggplot(df_counts) +
  geom_bar(aes(x = Value, y = Frequency, fill = Column), stat = "identity", position = position_dodge(width = 0.9)) +
  geom_segment(aes(x = as.numeric(Value) - 0.4, xend = as.numeric(Value) + 0.4,
                   y = CommonFreq, 
                   yend = CommonFreq, 
                   group = 1), color = "black") +
  labs(title = paste(title_name, ori_col_name, "and", ver1_col_name),
       x = "Value", y = "Frequency", fill = "Column") +
  theme_minimal() +
  scale_x_discrete(limits = all_values) + # Set all unique values for x-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Plot where the difference come from
df_diff <- df %>%
  filter(!!sym(ori_col_name) != !!sym(ver1_col_name)) %>%
  group_by_at(vars(ori_col_name, ver1_col_name)) %>%
  summarise(Frequency = n(), .groups = 'drop')

ggplot(df_diff,aes(x = !!sym(ori_col_name), y = !!sym(ver1_col_name), size = Frequency)) +
  geom_point(alpha = 0.7, color = "blue") +
  scale_x_discrete(limits = all_values) + # Set all unique values for x-axis
  scale_y_discrete(limits = all_values) + # Set all unique values for y-axis
  labs(title = paste(title_name, ori_col_name, "and", ver1_col_name),
       x = ori_col_name,
       y = ver1_col_name,
       size = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

## Study some problematic cases

col_name1 <- "R=100"
col_name2 <- "R=400"

ori_r100 = read.table(file = paste0(path_restest, 'original_R100', "/FinalTable.txt"))
ori = read.table(file = paste0(path_restest, 'original', "/FinalTable.txt"))
ver1_r100 = read.table(file = paste0(path_restest, 'ver1_R100', "/FinalTable.txt"))
ver1 = read.table(file = paste0(path_restest, 'ver1', "/FinalTable.txt"))

table(ver1_r100[suspect_case,"pred.y"])
table(ver1[suspect_case,"pred.y"])

suspect_case = which(ori_r100$pred.y == 22 & ori$pred.y ==1)
boxplot(Data_Res_Test[suspect_case, 5:9], outline = FALSE, ylim = c(-4,4), yaxt = 'n')
at_ticks <- seq(from = -4, to = 4, by = 2)  # Adjust the 'by' value as needed
axis(side = 2, at = at_ticks, las = 1)

boxplot(ori[which(ori$pred.y==1 & is.na(ori$Z.truth)), 5:9], outline = FALSE)
boxplot(Data_Res_Test[which(ori_r100$pred.y==3), 5:9], outline = FALSE)
boxplot(Data_Res_Test[which(ver1_r100$pred.y==1), 5:9], outline = FALSE)
boxplot(Data_Res_Test[which(ver1_r100$pred.y==3), 5:9], outline = FALSE)

config_c =35
con1 = ver1_r100$pred.y==config_c
con2 = ver1$pred.y==config_c
con3 = is.na(ver1_r100$Z.truth)
con4 = is.na(ver1$Z.truth)

df_plot = rbind(Data_Res_Test[which(con1&con3),],
                Data_Res_Test[which(con2&con4),]) %>%
  select(starts_with("T")) %>%
  mutate(version = c(rep("ver1100", length(which(con1&con3))), 
         rep("ver1400", length(which(con2&con4))))) 
df_plot$version = factor(df_plot$version)
df_long <- pivot_longer(df_plot[,-1], cols = starts_with("T"), names_to = "Variable", values_to = "Value")
df_long$Variable <- sub("Tvalue_", "", df_long$Variable)
df_long$Variable = factor(df_long$Variable, levels = name_six_diff[-1])
ggplot(df_long, aes(x = Variable, y = Value, color = version))+
  labs(title = paste0("Predicted configuration", config_c, " (not in the table)"))+
  geom_boxplot(outlier.shape = NA)+
  theme_minimal() +
  ylim(-5,5)+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text = element_text(size = 12))
