

predictiver_rule_ver4 <- function(significance_level,
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
  

  ### Remove insifnificant G-E
  remove_var = NA
  #' construct the logical table 
  #' 
  rm.ind = which(remove_var %in% list_name_test)
  length_data = nrow(Data_Res_Test)
  
  Z_trunc_final_code = construct_logical_table(prob, keep_config, remove_var, list_name_test,
                                               path_save = paste0(file_path_Results,version,"/"))
  
  List_names_tot = colnames(Z_trunc_final_code)
  if(length(rm.ind)!=0){
    List_names_final = List_names_tot[-rm.ind]
  }else{
    List_names_final = List_names_tot
  }
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
  Data_Res_Test  <- Data_Res_Test %>% 
    filter(abs(tGE)>1.96) 
  
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
    config_list_final[which(duplicated2(rbind(Z_trunc_final[,List_names_final], 
                                              Data_Res_coded[x, List_names_final])))[1]]
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
    
    
    while ((existing_pop_Learn!=17) & (existing_pop_Test !=17)){
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
            print("full")
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
              set_names(List_names_final) %>%
              bind_cols() %>% 
              as.data.frame()
            print("resampling")
            }
        } else {
          Nb_config_c <- R*0.2
          true_row <- intersect(which(Z_truth == config_true), c(1:length_data)[-trainIndex])
          if (length(true_row) > Nb_config_c){
            print("full")
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
              set_names(List_names_final) %>%
              bind_cols() %>% 
              as.data.frame()
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
  
  ### CONTINUE HERE 
  # apply the best rule to the real data  -----------------------------------
  
  RealData_x <- Data_Res_Test[,paste0("t",List_names_final)]
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

combine_rules <- function(version, 
                          B,
                          file_path_Results,
                          Data_Res_Test,
                          significance_level,
                          offset,
                          GE, 
                          number_pop){
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
  
  list_config = readRDS(paste0(file_path_Results,
                               version,
                               "/List_config.rds"))
  config_list_final = as.numeric(rownames(list_config))
  
  if(length(List_names_final)<6){
    RealData_x <- Data_Res_Test[,tail(names(Data_Res_Test), 5)]
  }else{
    RealData_x <- Data_Res_Test[,tail(names(Data_Res_Test), 6)]
  }
  colnames(RealData_x) <- List_names_final
  
  All_pred <- as.data.frame(matrix(NA,
                                   nrow = nrow(RealData_x),
                                   ncol = length(list_selected_rule)))
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
    RealData_predy_i <- predict(Model_i, newdata = RealData_x) 
    All_pred[,i] <- config_list_final[RealData_predy_i]
  }
  colnames(All_pred) <- paste0("iter", list_selected_rule)
  return(All_pred)
}

plot_similiar <- function(result, version, names_iter){
  names(result) <- c("iter1", "iter2")
  all_class = sort(unique(unlist(result, use.names = FALSE)))
  ratio_matrix <- matrix(0, nrow = length(all_class), ncol = length(all_class),
                         dimnames = list(all_class, all_class))
  
  for (val1 in 1:length(all_class)) {
    for (val2 in 1:length(all_class)) {
      freq_iter1 <- sum(result$iter1 == all_class[val1])
      freq_iter2 <- sum(result$iter2 == all_class[val2])
      common_freq <- sum(result$iter1 == all_class[val1] & result$iter2 == all_class[val2])
      ratio <- ifelse(max(freq_iter1, freq_iter2) > 0, common_freq / max(freq_iter1, freq_iter2), 0)
      ratio_matrix[val1, val2] <- ratio
    }
  }
  
  melted_matrix <- reshape2::melt(ratio_matrix, varnames = c("Value in iter1", "Value in iter2"))
  
  # Step 4: Plot the heatmap
  p <- p <- ggplot(melted_matrix, aes(x = `Value in iter1`, y = `Value in iter2`, fill = value)) +
    geom_tile(color = "grey", size = 0.1) +  # Add tile borders for distinction
    scale_fill_gradient(low = "white", high = "darkgreen") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4.5),
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size = 4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.2, "cm"),
          legend.title = element_blank(),
          plot.tag = element_text(size = 6),
          plot.subtitle = element_text(size = 6),
          legend.box.spacing = unit(0, "pt"),
          strip.text = element_text(size = 4),
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.background = element_rect(fill = "white")) +
    xlab(names_iter[1]) + 
    ylab(names_iter[2]) +
    scale_x_discrete(limits = all_class) +
    scale_y_discrete(limits = all_class)
  
  jpeg(paste0(path_results,"figure/", version, names_iter[1], names_iter[2], ".jpeg"),
       width = 2000, height = 1700, res = 600) # chan
  print(p)
  # print(ml1)
  dev.off() 
}

