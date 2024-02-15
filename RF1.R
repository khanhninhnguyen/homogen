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
B = 2
offset=0
GE=0
number_pop = 3
R = 100
prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
          0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
keep_config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
version = "version_six_tests"

remove_var = "G-E"
list_name_test = c("G-E", "G-G'", "G-E'", "E-E'", "G'-E'","G'-E")

library(caret)
source(file = paste0(path_code, "support_RF1.R"))
source(file = paste0(path_code, "test_predictive_rule.R"))

path_restest <- paste0(path_results,"attribution/predictive_rule/")
file_path_Results=paste0(path_results,'attribution/predictive_rule/')


#' read test results 
List_main = read.table(file = paste0(path_restest,"list_selected_nmin200_10nearby.txt"), 
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

Data_Res_Test_fillNA <- Data_Res_Test %>%
  arrange(main, brp) %>%
  group_by(main, brp) %>%
  fill(Tvalue_GPS_ERA, .direction = "downup") %>%
  ungroup()
#' remove the err cases
#' 
find_bug <- function(main_beg_new, main_end_new, nearby_beg_new, nearby_end_new){
  cond = as.integer((main_beg_new > nearby_beg_new) & 
                      (main_end_new > nearby_end_new))
  return(cond)
}
fix_case <- List_main %>% 
  rowwise() %>%
  mutate(Fix = find_bug(main_beg_new,
                        main_end_new,
                        nearby_beg_new,
                        nearby_end_new))

# Data_Res_Test <- Data_Res_Test[which(fix_case$Fix == 0),]
# rownames(Data_Res_Test) <- NULL

# List_main <- List_main[which(fix_case$Fix == 0),]
# rownames(List_main) <- NULL

# Data_Res_Test_fillNA <- Data_Res_Test_fillNA[which(fix_case$Fix == 0),]
Data_Res_Test <- Data_Res_Test_fillNA[which(fix_case$Fix == 0),]
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
#                   Data_Res_Test, 
#                   path_restest,
#                   version = "original")
# write.table(a, file = paste0(path_restest, 'original', "/FinalTable.txt"), sep = '\t', quote = FALSE)

a1 = predictiver_rule_ver2(significance_level, B, 
                              offset, 
                              GE,
                              number_pop,
                              R,
                              prob, 
                              keep_config,
                              remove_var,
                              list_name_test,
                              Data_Res_Test, 
                              # Data_Res_Test_fillNA,
                              path_restest,
                              version = "ver1_R100")
write.table(a1, file = paste0(path_restest, 'ver1_R100', "/FinalTable.txt"), sep = '\t', quote = FALSE)



