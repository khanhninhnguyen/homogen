#' test to find the best predictive rule 
#' ##### function ---------


predictiver_rule_ver2 <- function(significance_level,
                                  B, 
                                  offset,
                                  GE,
                                  number_pop,
                                  R,
                                  prob, 
                                  keep_config,
                                  remove_var,
                                  list_name_test,
                                  Data_Res_Test, 
                                  path_restest,
                                  version){
  
  
  #' construct the logical table 
  #' 
  rm.ind = which(remove_var %in% list_name_test)
  length_data = nrow(Data_Res_Test)
  
  Z_trunc_final_code = construct_logical_table(prob, keep_config, remove_var, list_name_test,
                                               path_save = path_restest)
  List_names_tot = colnames(Z_trunc_final_code)
  List_names_final = List_names_tot[-rm.ind]
  print(List_names_final)
  config_list_final = as.numeric(rownames(Z_trunc_final_code))
  
  # Reverse the transformations
  Z_trunc_final <- Z_trunc_final_code  # Start with a copy
  Z_trunc_final[Z_trunc_final_code == 3] <- 1
  Z_trunc_final[Z_trunc_final_code == 2] <- -1
  Z_trunc_final[Z_trunc_final_code == 1] <- 0
  
  Nbconfig <- nrow(Z_trunc_final_code)
  NbSim <- R * Nbconfig
  
  colnames(Data_Res_Test)[4:9] <- paste0("t", List_names_tot) 
  
  #' check the configuration based on test only 
  t_cols <- colnames(Data_Res_Test)[4:9]
  
  # Perform the operations
  Data_Res_coded <- Data_Res_Test %>%
    mutate(across(all_of(t_cols),
                  ~ {
                    # Step 1: Convert t-value to p-value
                    p_value <- 2 * pnorm(-abs(.x))
                    # Step 2: Convert p-value to 0 or 1 based on threshold
                    p_converted <- if_else(p_value < 0.05, 1, 0)
                    # Step 3: Multiply by sign of t-value
                    result <- sign(.x) * p_converted
                    return(result)
                  },
                  .names = "code_{.col}")) %>%
  select(starts_with("code_"))  %>%
  rename_with(.fn = ~ sub("code_t", "", .x), .cols = starts_with("code_"))
  
  # find the truth based on test 
  Z_truth <- sapply(c(1:nrow(Data_Res_coded)), function(x) {
    config_list_final[which(duplicated2(rbind(Z_trunc_final[,-rm.ind], 
                            Data_Res_coded[x, -rm.ind])))[1]]
  })
  
  results_list <- lapply(setNames(List_names_final, List_names_final), function(name) {
    Keep.Data(name, significance.level = significance_level, Data.Res.Test = Data_Res_Test)
  })
  
  for (name in names(results_list)) {
    assign(paste0("Data_", name), results_list[[name]])
  }
  
  # run the bootstrap 
  
  set.seed(1)
  
  error_RF <- rep(NA, B)
  
  for (b in 1:B){
    
    print(b)
    
    ######
    # Vfold 
    existing_pop_Learn <- 0
    existing_pop_Test <- 0
    
    
    while ((existing_pop_Learn!=15) & (existing_pop_Test !=15)){
      trainIndex <- createDataPartition(1:length_data, p = 0.8, list = FALSE)
      
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
      calculate_unique_pop_sum <- function(phase, dataset_names) {
        sum(sapply(dataset_names, function(name) {
          data_with_name <- paste("Data", name, phase, sep = "_")
          dataset <- get(data_with_name)
          length(unique(dataset$pop))
        }))
      }
      
      existing_pop_Learn <- calculate_unique_pop_sum("Learn", List_names_final)
      existing_pop_Test <- calculate_unique_pop_sum("Test", List_names_final)
    }
    
    #####
    # Construction of the pop of t.values for each series Learn and Test
    for (i in 1:length(List_names_final)){
      eval(parse(text=paste0("Pop_",List_names_final[i],"_Learn=Pop.create(Data_",List_names_final[i],"_Learn)")))
      eval(parse(text=paste0("Pop_",List_names_final[i],"_Test=Pop.create(Data_",List_names_final[i],"_Test)")))
    }
    
    print("err")
    #####
    # Bootstrap: ccontruction de DataLearn et DataTest respecting the proportion of configuration use Boot1, 
    # Boot() for the same number of configuration
    DataLearn <- c()
    DataTest <- c()
    
    Boot <- function(type_dataset, Z_trunc_code, NbSim, List_names_final){
      config=c()
      Data.res <- c()
      Nbconfig <- nrow(Z_trunc_code)
      
      for (c in 1:Nbconfig){
        print(c)
        config_code <- Z_trunc_code[c,]
        config_true = as.numeric(rownames(config_code))

        # N <- eval(parse(text=paste0("NbSim",type_dataset)))
        if (type_dataset=="Learn"){
          Nb_config_c <- R*0.8
          true_row <- intersect(which(Z_truth == config_true), trainIndex)
          if (length(true_row) > Nb_config_c){
            ind_sel <- sample(true_row, Nb_config_c)
            data_c <- Data_Res_Test[ind_sel, paste0("t", List_names_final)]
            colnames(data_c) <- sub("^t", "", colnames(data_c))
          }else{
            data_c <- purrr::map(List_names_final,~{
              Pop <- c()
              code_names_series <- config_code[which(colnames(Z_trunc_code) %in% .x)]
              eval(parse(text = paste0("Pop=Pop_",.x,"_",type_dataset,"[[",code_names_series,"]]")))
              res.t = sample(Pop, Nb_config_c, replace = FALSE) # maybe change here 
              return(res.t)
            }) %>% 
              bind_cols() %>% 
              as.data.frame()
            colnames(data_c) = List_names_final
          }
        } else {
          Nb_config_c <- R*0.2
          true_row <- intersect(which(Z_truth == config_true), c(1:length_data)[-trainIndex])
          if (length(true_row) > Nb_config_c){
            ind_sel <- sample(true_row, Nb_config_c)
            data_c <- Data_Res_Test[ind_sel, paste0("t", List_names_final)]
            colnames(data_c) <- sub("^t", "", colnames(data_c))
          }else{
            data_c <- purrr::map(List_names_final,~{
              Pop <- c()
              code_names_series <- config_code[which(colnames(Z_trunc_code) %in% .x)]
              eval(parse(text = paste0("Pop=Pop_",.x,"_",type_dataset,"[[",code_names_series,"]]")))
              res.t = sample(Pop, Nb_config_c, replace = FALSE) # maybe change here 
              return(res.t)
            }) %>% 
              bind_cols() %>% 
              as.data.frame()
            colnames(data_c) = List_names_final
          }
        }
        
        config <- c(config,rep(c,Nb_config_c))
        Data.res = rbind(Data.res,data_c)
      }
      
      Data.res <- Data.res %>% mutate(config=config)
      Data.res$config <- as.factor(Data.res$config)
      
      return(Data.res)
    }
    
    
    DataLearn <- Boot("Learn",Z_trunc_final_code, NbSim, List_names_final)
    DataTest <- Boot("Test",Z_trunc_final_code, NbSim, List_names_final)
    
    saveRDS(DataLearn, file = paste0(
      file_path_Results, version,
      "/DataLearn_", 
      b,
      significance_level, 
      offset, 
      GE, 
      number_pop,
      ".rds"))
    saveRDS(DataTest, file = paste0(
      file_path_Results, version,
      "/DataTest_",
      b,
      significance_level,
      offset, 
      GE,
      number_pop,
      ".rds"))
    
    Res.pred1 <- PredRule_RF(DataLearn, DataTest, b, Nbconfig)
    saveRDS(Res.pred1, file = paste0(
      file_path_Results, version,
      "/Res.pred_",
      b,
      significance_level, 
      offset, 
      GE, 
      number_pop,
      ".rds"))
    error_RF[b] <- Res.pred1$err.tot 
  }
  print(error_RF)
  # read the best predictive rule
  
  FinalPred <- readRDS(paste0(file_path_Results, version,
                              "/Res.pred_",
                              which.min(error_RF),
                              significance_level, 
                              offset, 
                              GE, 
                              number_pop,
                              ".rds"))
  FinalModel <- FinalPred$modrf
  
  # apply the best rule to the real data  -----------------------------------

  RealData_x <- Data_Res_Test[,colnames(Data_Res_Test) %in% c("tGGp","tGEp","tEEp", "tGpEp","tGpE")]
  colnames(RealData_x) <- List_names_final
  RealData_predy <- predict(FinalModel, newdata = RealData_x) 
  Final_pred_y <- config_list_final[RealData_predy]
  
  # synthesize results 
  colnames(Data_Res_coded) <- paste0("code.",colnames(Data_Res_coded))
  FinalTable <- Data_Res_Test %>%
    cbind(Data_Res_coded) %>%
    rename_with(~str_replace(.x, "^t", ""), starts_with("t")) %>%
    mutate(Z.truth = Z_truth,
           pred.y = Final_pred_y)
  
  return(FinalTable)
  
}


predictiver_rule_original <- function(significance_level,
                                      B, 
                                      offset,
                                      GE,
                                      number_pop,
                                      R,
                                      prob, 
                                      keep_config,
                                      remove_var,
                                      list_name_test,
                                      Data_Res_Test, 
                                      path_restest,
                                      version){
  
  
  #' construct the logical table 
  #' 
  Z_trunc_final_code = construct_logical_table(prob, keep_config, remove_var, list_name_test,
                                               path_save = path_restest)
  List_names_tot = colnames(Z_trunc_final_code)
  Nbconfig <- nrow(Z_trunc_final_code)
  NbSim <- R * Nbconfig
  
  colnames(Data_Res_Test)[4:9] <- paste0("t", List_names_tot) 
  
  # Separate into learning and test dataset 
  
  rm.ind = which(remove_var %in% list_name_test)
  List_names_final = List_names_tot[-rm.ind]
  print(List_names_final)
  
  results_list <- lapply(setNames(List_names_final, List_names_final), function(name) {
    Keep.Data(name, significance.level = significance_level, Data.Res.Test = Data_Res_Test)
  })
  
  for (name in names(results_list)) {
    assign(paste0("Data_", name), results_list[[name]])
  }
  
  # run the bootstrap 
  
  set.seed(1)
  
  error_RF <- rep(NA, B)
  
  for (b in 1:B){
    
    print(b)
    
    ######
    # Vfold 
    existing_pop_Learn <- 0
    existing_pop_Test <- 0
    
    
    while ((existing_pop_Learn!=15) & (existing_pop_Test !=15)){
      trainIndex <- createDataPartition(1:nrow(Data_Res_Test), p=0.8, list = FALSE)
      
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
      calculate_unique_pop_sum <- function(phase, dataset_names) {
        sum(sapply(dataset_names, function(name) {
          data_with_name <- paste("Data", name, phase, sep = "_")
          dataset <- get(data_with_name)
          length(unique(dataset$pop))
        }))
      }
      
      existing_pop_Learn <- calculate_unique_pop_sum("Learn", List_names_final)
      existing_pop_Test <- calculate_unique_pop_sum("Test", List_names_final)
    }
    
    #####
    # Construction of the pop of t.values for each series Learn and Test
    for (i in 1:length(List_names_final)){
      eval(parse(text=paste0("Pop_",List_names_final[i],"_Learn=Pop.create(Data_",List_names_final[i],"_Learn)")))
      eval(parse(text=paste0("Pop_",List_names_final[i],"_Test=Pop.create(Data_",List_names_final[i],"_Test)")))
    }
    
    print("err")
    #####
    # Bootstrap: ccontruction de DataLearn et DataTest respecting the proportion of configuration use Boot1, 
    # Boot() for the same number of configuration
    DataLearn <- c()
    DataTest <- c()
    
    Boot <- function(type.dataset, Z.trunc.code, NbSim, List.names.final){
      config=c()
      Data.res <- c()
      Nbconfig <- nrow(Z.trunc.code)
      
      for (c in 1:Nbconfig){
        print(c)
        config.code <- Z.trunc.code[c,]
        # N <- eval(parse(text=paste0("NbSim",type.dataset)))
        if (type.dataset=="Learn"){
          Nb.config.c <- R*0.8
        } else {Nb.config.c <- R*0.2}
        config <- c(config,rep(c,Nb.config.c))
        
        data.c <- purrr::map(List.names.final,~{
          Pop <- c()
          code.names.series <- config.code[which(colnames(Z.trunc.code) %in% .x)]
          eval(parse(text=paste0("Pop=Pop_",.x,"_",type.dataset,"[[",code.names.series,"]]")))
          res.t=sample(Pop, Nb.config.c, replace = FALSE)
          return(res.t)
        }) %>% bind_cols() %>% as.data.frame()
        colnames(data.c)=List.names.final
        Data.res = rbind(Data.res,data.c)
      }
      
      Data.res <-Data.res %>% mutate(config=config)
      Data.res$config <- as.factor(Data.res$config)
      
      return(Data.res)
    }
    
    
    DataLearn <- Boot("Learn",Z_trunc_final_code, NbSim, List_names_final)
    DataTest <- Boot("Test",Z_trunc_final_code, NbSim, List_names_final)
    
    saveRDS(DataLearn, file = paste0(
      file_path_Results, version,
      "/DataLearn_", 
      b,
      significance_level, 
      offset, 
      GE, 
      number_pop,
      ".rds"))
    saveRDS(DataTest, file = paste0(
      file_path_Results, version,
      "/DataTest_",
      b,
      significance_level,
      offset, 
      GE,
      number_pop,
      ".rds"))
    
    Res.pred1 <- PredRule_RF(DataLearn, DataTest, b, Nbconfig)
    saveRDS(Res.pred1, file = paste0(
      file_path_Results, version,
      "/Res.pred_",
      b,
      significance_level, 
      offset, 
      GE, 
      number_pop,
      ".rds"))
    error_RF[b] <- Res.pred1$err.tot 
  }
  print(error_RF)
  # read the best predictive rule
  
  FinalPred <- readRDS(paste0(file_path_Results, version,
                              "/Res.pred_",
                              which.min(error_RF),
                              significance_level, 
                              offset, 
                              GE, 
                              number_pop,
                              ".rds"))
  FinalModel <- FinalPred$modrf
  
  # apply the best rule to the real data  -----------------------------------
  
  # check cases in the table
  p.values.i=c()
  truth.vec.i=c()
  Z.truth.i <- c()
  
  config_list_final = as.numeric(rownames(Z_trunc_final_code))
  Z_trunc_final <- Z_trunc_final_code  # Start with a copy
  
  # Reverse the transformations
  Z_trunc_final[Z_trunc_final_code == 3] <- 1
  Z_trunc_final[Z_trunc_final_code == 2] <- -1
  Z_trunc_final[Z_trunc_final_code == 1] <- 0
  for (i in 1:nrow(Data_Res_Test)){
    a <- c()
    p.values.i <- 2*pnorm(-abs(as.numeric(Data_Res_Test[i,paste0("t", List_names_final)])))
    a <- ifelse(p.values.i < significance_level, 1, 0)*sign(as.numeric(Data_Res_Test[i,paste0("t", List_names_final)]))
    truth.vec.i <- rbind(truth.vec.i,a)
    Z.truth.i <- c(Z.truth.i,config_list_final[which(duplicated2(rbind(Z_trunc_final[,-rm.ind],a)))[1]])
  }
  
  pred.truth <- as.data.frame(cbind(truth.vec.i,Z.truth.i))
  colnames(pred.truth) <- c("code.GGp", "code.GEp" , "code.EEp", "code.GpEp","code.GpE","Z.truth")
  
  Data_Res_Test <- cbind(Data_Res_Test,pred.truth)
  RealData.x <- Data_Res_Test[,colnames(Data_Res_Test) %in% c("tGGp","tGEp","tEEp", "tGpEp","tGpE")]
  colnames(RealData.x) <- List_names_final
  
  RealData.predy <- predict(FinalModel,newdata=RealData.x) 
  FinalTable <- cbind(Data_Res_Test,config_list_final[RealData.predy])
  colnames(FinalTable)[which(colnames(FinalTable)=="config_list_final[RealData.predy]")] <- "pred.y"
  
  return(FinalTable)
  
}

