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
B = 50
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

# remove_var = "G-E"
list_name_test = c("G-E", "G-G'", "G-E'", "E-E'", "G'-E'","G'-E")

library(caret)
source(file = paste0(path_code, "support_RF1.R"))
source(file = paste0(path_code, "test_predictive_rule.R"))

path_restest <- paste0(path_results,"attribution/FGLS_comb/")
file_path_Results=paste0(path_results,'attribution/predictive_rule/')


# read test results  -------------------
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


# run the training of predictive rule -------------------------------------
version = "ver1/400/"

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
                              version = "ver4/")
write.table(a1, file = paste0(file_path_Results, "ver4", "/FinalTable.txt"), sep = '\t', quote = FALSE)


# give the final models to predict ------------------------

version = "ver4"
error_list = rep(NA,B)

for (b in c(1:B)) {
  FinalPred <- readRDS(paste0(file_path_Results,
                              version,
                              "/Res.pred_",
                              b,
                              significance_level, 
                              offset, 
                              GE, 
                              number_pop,
                              ".rds"))
  error_list[b] = FinalPred$err.tot
}
List_names_final = FinalPred$modrf$coefnames
list_selected_rule = which(error_list == min(error_list))
list_selected_rule

all_model <- list()
set.seed(1)
for (i in 1:length(list_selected_rule)) {
  FinalPred_i <- readRDS(paste0(file_path_Results,
                                version,
                                "/Res.pred_",
                                list_selected_rule[i],
                                significance_level, 
                                offset, 
                                GE, 
                                number_pop,
                                ".rds"))
  Model_i <- FinalPred_i$modrf
  all_model[[paste0("model",i)]] <- Model_i
}

ind <- rep(NA, 24)
Res <- rep(NA, 24)

set.seed(1)

for (j in c(1:24)) {
  Modi = all_model[[paste0("model",list_selected_rule[j])]] 
  ind[j] = predict(Modi, newdata = X)
  Res[j] = config_list_final[ind[j]]
}
Res
table(Res)
table(ind)

save(all_model, file = paste0(path_results, "RF_models.RData"))

# list_config = readRDS(paste0(file_path_Results,
#                              version,
#                              "/List_config.rds"))
# config_list_final = as.numeric(rownames(list_config))
# 
# if(length(List_names_final)<6){
#   RealData_x <- Data_Res_Test[,tail(names(Data_Res_Test), 5)]
# }else{
#   RealData_x <- Data_Res_Test[,tail(names(Data_Res_Test), 6)]
# }
# colnames(RealData_x) <- List_names_final
# 
# All_pred <- as.data.frame(matrix(NA,
#                                  nrow = nrow(RealData_x),
#                                  ncol = length(list_selected_rule)))
# for (i in 1:length(list_selected_rule)) {
#   FinalPred_i <- readRDS(paste0(file_path_Results,
#                                 version,
#                                 "/Res.pred_",
#                                 list_selected_rule[i],
#                                 significance_level, 
#                                 offset, 
#                                 GE, 
#                                 number_pop,
#                                 ".rds"))
#   Model_i <- FinalPred_i$modrf
#   RealData_predy_i <- predict(Model_i, newdata = RealData_x[1,]) 
#   All_pred[,i] <- config_list_final[RealData_predy_i]
# }



# apply the final models to test results  ---------------------------------

pred_rules = get(load( file = paste0(path_results, "RF_models.RData")))

X <- Data_Res_Test[,4:9] %>%
  setNames( c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE"))
  

for (i in c(1:length(pred_rules))) {
  
  Model_i = pred_rules[[paste0("model",i)]]
  predict_i = predict(Model_i, newdata = X)
  
  Data_Res_Test[paste("iter", i)] = predict_i
  
}









