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
R = 10
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

predictiver_rule <- function(significance_level, 
                             B = 20, 
                             offset, 
                             GE,
                             number_pop,
                             R,
                             prob, 
                             keep_config,
                             remove_var,
                             list_name_test, 
                             Data_Res_Test, 
                             path_restest){
  
  
  #' construct the logical table 
  #' 
  Z_trunc_final_code = construct_logical_table(prob, keep_config, remove_var, list_name_test,
                                               path_save = path_restest)
  List_names_tot = colnames(Z_trunc_final_code)
  Nbconfig <- nrow(Z_trunc_final_code)
  NbSim <- R * Nbconfig
  
  colnames(Data_Res_Test)[4:9] <- paste0("t", List_names_tot) 
 
  # Separate into learning and test dataset 
  
  List_names_final = List_names_tot[-which(remove_var %in% list_name_test)]
  for (i in  1:length(List_names_final)){
    eval(parse(text=paste0("Data_",
                           List_names_final[i],
                           "=Keep.Data(List_names_final[i], significance.level = ", 
                           significance_level, 
                           ", Data.Res.Test = ",
                           Data_Res_Test,
                           ")")))  
  }
  
  # run the bootstrap 
  
  set.seed(1)
  
  error_RF <- rep(NA, B)
  
  for (b in 1:B){
    
    ######
    # Vfold 
    existing_pop_Learn <- 0
    existing_pop_Test <- 0
    
    
    while ((existing_pop_Learn!=15) & (existing_pop_Test !=15)){
      trainIndex <- createDataPartition(1:nrow(Data.GGp) , p=0.8, list = FALSE)
      
      # Loop through each dataset name
      for (name in List_names_final) {
        # Construct the variable names for learning and testing sets
        learn_name <- paste("Data", name, "Learn", sep = "_")
        test_name <- paste("Data", name, "Test", sep = "_")
        
        # Slice the data for learning and testing sets
        assign(learn_name, slice(get(paste("Data", name, sep = "_")), trainIndex))
        assign(test_name, slice(get(paste("Data", name, sep = "_")), -trainIndex))
      }
      
      # check if each dataset contain enough 3 population
      existing_pop_Learn <- calculate_unique_pop_sum("Learn", List_names_final)
      existing_pop_Test <- calculate_unique_pop_sum("Test", List_names_final)
    }
    
    #####
    # Construction of the pop of t.values for each series Learn and Test
    for (i in 1:length(List_names_final)){
      eval(parse(text=paste0("Pop.",List_names_final[i],".Learn=Pop.create(Data.",List_names_final[i],".Learn)")))
      eval(parse(text=paste0("Pop.",List_names_final[i],".Test=Pop.create(Data.",List_names_final[i],".Test)")))
    }
    
    #####
    # Bootstrap: ccontruction de DataLearn et DataTest respecting the proportion of configuration use Boot1, 
    # Boot() for the same number of configuration
    DataLearn <- c()
    DataTest <- c()
    
    DataLearn <- Boot("Learn",Z_trunc_final_code, NbSim, List_names_final)
    DataTest <- Boot("Test",Z_trunc_final_code, NbSim, List_names_final)
    
    saveRDS(DataLearn, file = paste0(
      file_path_Results,
      "DataLearn_", 
      b,
      significance_level, 
      offset, 
      GE, 
      number_pop,
      ".rds"))
    saveRDS(DataTest, file = paste0(
      file_path_Results,
      "DataTest_",
      b,
      significance_level,
      offset, 
      GE,
      number_pop,
      ".rds"))
    
    Res.pred1 <- PredRule_RF(DataLearn,DataTest,b,Nbconfig)
    saveRDS(Res.pred1, file = paste0(
      file_path_Results,
      "Res.pred_",
      b,
      significance_level, 
      offset, 
      GE, 
      number_pop,
      ".rds"))
    error_RF[b] <- Res.pred1$err.tot 
    print(b)
  }
  
  # read the best predictive rule
  
  FinalPred <- readRDS(paste0(file_path_Results,
                              "modrf_b",
                              b = which.min(error_RF),
                              significance_level, 
                              offset, GE, 
                              number_pop,
                              ".rds"))
  
}

