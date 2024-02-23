#' Train classifier RANDOM FOREST
#' input:
#' * test results from FGLS 
#' output:
#' list of changepoints attributed 
#' Final predictive rule program 
#' 
#' 
####### Test result - remove when finish ##################
  
significance_level = 0.05
B = 20
offset=0
GE=0
number_pop = 3
R = 400
prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
          0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
keep_config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
version = "ver1/400/"

remove_var = "G-E"
list_name_test = c("G-E", "G-G'", "G-E'", "E-E'", "G'-E'","G'-E")

library(caret)
source(file = paste0(path_code, "support_RF1.R"))
source(file = paste0(path_code, "test_predictive_rule.R"))

path_restest <- paste0(path_results,"attribution/FGLS_comb/")
file_path_Results=paste0(path_results,'attribution/predictive_rule/')


#' read test results 
List_main = read.table(file = paste0(path_results,"list_selected_nmin200_10nearby.txt"), 
                       header = TRUE, 
                       stringsAsFactors = FALSE) 

name.version ="FGLS_jump_tvalue.txt"
name.results <- paste0(path_restest, name.version) # name of test result file
Data_Res_Test0 <- read.table(name.results,
                             header = TRUE,
                             stringsAsFactors = FALSE) 

Data_Res_Test <- cbind(List_main[,c("main", "brp", "nearby")],
                       Data_Res_Test0[,7:12]) %>%
  mutate(brp = as.Date(List_main$brp, format="%Y-%m-%d")) 

Data_Res_Test <- Data_Res_Test %>%
  arrange(main, brp) %>%
  group_by(main, brp) %>%
  fill(Tvalue_GPS_ERA, .direction = "downup") %>%
  ungroup() %>%
  filter(abs(Tvalue_GPS_ERA)>1.96)

#' remove the err cases
#' 
# find_bug <- function(main_beg_new, main_end_new, nearby_beg_new, nearby_end_new){
#   cond = as.integer((main_beg_new > nearby_beg_new) & 
#                       (main_end_new > nearby_end_new))
#   return(cond)
# }
# fix_case <- List_main %>% 
#   rowwise() %>%
#   mutate(Fix = find_bug(main_beg_new,
#                         main_end_new,
#                         nearby_beg_new,
#                         nearby_end_new))

# Data_Res_Test <- Data_Res_Test[which(fix_case$Fix == 0),]
# rownames(Data_Res_Test) <- NULL

# List_main <- List_main[which(fix_case$Fix == 0),]
# rownames(List_main) <- NULL

# Data_Res_Test_fillNA <- Data_Res_Test_fillNA[which(fix_case$Fix == 0),]
# 
# a = predictiver_rule_original(significance_level, B,
#                   offset,
#                   GE,
#                   number_pop,
#                   R,
#                   prob,
#                   keep_config,
#                   remove_var,
#                   list_name_test,
#                   Data_Res_Test_fillNA,
#                   path_restest,
#                   version = "original/R400/")
# write.table(a, file = paste0(file_path_Results, 'original/R400/', "/FinalTable.txt"), sep = '\t', quote = FALSE)

a = predictiver_rule_ver2(significance_level, 
                          B,
                          offset,
                          GE,
                          number_pop,
                          R,
                          prob,
                          keep_config,
                          remove_var = "G-E",
                          list_name_test,
                          Data_Res_Test_fillNA,
                          path_restest,
                          version = "ver1/R400/")
write.table(a, file = paste0(file_path_Results, 'ver1/R400/', "/FinalTable.txt"), sep = '\t', quote = FALSE)

a1 = predictiver_rule_ver4(significance_level, B, 
                              offset, 
                              GE,
                              number_pop,
                              R,
                              prob, 
                              keep_config,
                              remove_var = NA,
                              list_name_test,
                              Data_Res_Test, 
                              path_restest,
                              version = "ver4b/")
write.table(a1, file = paste0(file_path_Results, "ver4b", "/FinalTable.txt"), sep = '\t', quote = FALSE)

# check the similarity between different iteration ------------------------
all_pred = combine_rules(version = "ver4", 
                        B,
                        file_path_Results,
                        Data_Res_Test,
                        significance_level,
                        offset,
                        GE, 
                        number_pop)
Final_table = read.table(file = paste0(file_path_Results, version = "ver4", "/FinalTable.txt"))
z_truth = Final_table$Z.truth
# all_pred$truth = z_truth
# all_pred$iden = apply(all_pred, 1, function(row) {
#   if (is.na(row['truth'])) {
#     return(NA)  # Return NA if the 'truth' value is NA
#   }
#   sum(row[-length(row)] == row['truth'], na.rm = TRUE)
# })
# all_pred$iden = all_pred$iden /7
# table(all_pred$iden)

all_pred$mean = all_res$voted_value
for (j in c(2:ncol(all_pred))) {
  plot_similiar(result = all_pred[which(is.na(z_truth)),c(1,j)], version = "ver4", names_iter = colnames(all_pred)[c(1,j)])
}

#' consistency between nearby stations 
all_res = cbind(Data_Res_Test[,c(1:3)], all_pred) 

check_same_value <- function(row) {
  as.numeric(length(unique(row)) == 1)
}

find_single_mode <- function(row) {
  freq <- table(row)
    return(names(freq)[which.max(freq)])
}

# Apply the function across the iter columns
all_res$same_value <- apply(all_res[,4:10], 1, check_same_value)
all_res$voted_value <- apply(all_res[,4:10], 1, find_single_mode)
all_res$truth = z_truth
all_res$dd = List_main$dd


