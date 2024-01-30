#' Train classifier RANDOM FOREST
#' input:
#' * test results from FGLS 
#' output:
#' list of changepoints attributed 
#' Final predictive rule program 
#' 
#' 
####### Test result - remove when finish ##################
significance.level = 0.05
B = 20
offset=0
GE=0
number.pop = 3
R = 100
prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
          0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
keep_config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)

remove_var = "G-E"
list_name_test = c("G-E", "G-G'", "G-E'", "E-E'", "G'-E'","G'-E")

library(caret)
source(file = paste0(path_code, "support_RF1.R"))
path_restest <- paste0(path_results,"attribution/predictive_rule/")
file_path_Results=paste0(path_results,'attribution/predictive_rule/')

#' construct the logical table 
#' 
Z_trunc_final_code = construct_logical_table(prob, keep_config, remove_var, list_name_test,
                                             path_save = path_restest)
List_names_tot = colnames(Z_trunc_final_code)
Nbconfig <- nrow(Z_trunc_final_code)
NbSim <- R * Nbconfig

#' read test results 
List_main = read.table(file = paste0(path_restest,"list_selected_nmin200_10nearby.txt"), 
                       header = TRUE, 
                       stringsAsFactors = FALSE) 

name.version ="FGLS_jump_tvalue.txt"
name.results <- paste0(path_restest, name.version) # name of test result file
Data.Res.Test0 <- read.table(name.results,
                             header = TRUE,
                             stringsAsFactors = FALSE) 

Data.Res.Test <- cbind(List_main[,c("main", "brp", "nearby")],
                       Data.Res.Test0[,7:12]) %>%
  mutate(brp = as.Date(List_main$brp, format="%Y-%m-%d")) 

colnames(Data.Res.Test)[4:9] <- paste0("t", List_names_tot) 

# Separate into learning and test dataset 

List.names.final = List.names.tot[-rm.ind]
for (i in  1:length(List.names.final)){
  eval(parse(text=paste0("Data.",List.names.final[i],"=Keep.Data(List.names.final[i])")))  
}

# run the boostrap 

set.seed(1)

error.test <- rep(NA, B)





