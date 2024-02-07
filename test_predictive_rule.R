#' test to find the best predictive rule 
#' ##### function ---------
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
    
    print(Data_GGp_Learn)
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
          res.t=sample(Pop, Nb.config.c, replace = TRUE)
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
    
    Res.pred1 <- PredRule_RF(DataLearn, DataTest, b, Nbconfig)
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
  }
  
  # read the best predictive rule
  
  FinalPred <- readRDS(paste0(file_path_Results,
                              "Res.pred_",
                              b,
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
                                      path_restest){
  
  
  #' construct the logical table 
  #' 
  Z_trunc_final_code = construct_logical_table(prob, keep_config, remove_var, list_name_test,
                                               path_save = path_restest)
  List_names_tot = colnames(Z_trunc_final_code)
  Nbconfig <- nrow(Z_trunc_final_code)
  NbSim <- R * Nbconfig
  
  colnames(Data_Res_Test)[4:9] <- paste0("t", List_names_tot) 
  
  #' check the configuration based on test only 
  Data_Res_Test_coded <- Data_Res_Test %>%
    # Step 1: Calculate p-values for t-value columns
    mutate(across(starts_with("t"), 
                  ~ 2 * pnorm(-abs(.x)),
                  .names = "p_{.col}")) %>%
    # Step 2: Categorize p-values
    mutate(across(starts_with("p_"), 
                  ~ case_when(
                    . < 0.05 ~ 1,
                    TRUE ~ 0),
                  .names = "categorized_{.col}")) %>%
  ###### STOP here  
  
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
    
    print(Data_GGp_Learn)
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
          res.t=sample(Pop, Nb.config.c, replace = TRUE)
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
    
    Res.pred1 <- PredRule_RF(DataLearn, DataTest, b, Nbconfig)
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
  }
  
  # read the best predictive rule
  
  FinalPred <- readRDS(paste0(file_path_Results,
                              "Res.pred_",
                              b,
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

predictiver_rule_ver3 <- function(significance_level,
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
    
    print(Data_GGp_Learn)
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
          res.t=sample(Pop, Nb.config.c, replace = TRUE)
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
    
    Res.pred1 <- PredRule_RF(DataLearn, DataTest, b, Nbconfig)
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
  }
  
  # read the best predictive rule
  
  FinalPred <- readRDS(paste0(file_path_Results,
                              "Res.pred_",
                              b,
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

