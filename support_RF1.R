#' functions used in RF1
#' Construct the logical table 
#' @param prob probability of all configuration
#' @param keep_config index of configuration that want to be kept

prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
          0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)

keep_config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)

remove_var = "G-E"


construct_logical_table <- function(prob, keep_config, remove_var, list_name_test){
  G = c(rep(1,9), rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9))
  E = c(rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9),rep(1,9))
  Gp = rep(c(rep(0,3),rep(1,3),rep(-1,3)),6)
  Ep = rep(c(0,1,-1),18)
  
  Y = data.frame(G=G, E=E, Gp=Gp, Ep=Ep)
  Z = data.frame(GE=Y$G-Y$E, GGp=Y$G-Y$Gp, GEp=Y$G-Y$Ep,
                 EEp=Y$E-Y$Ep, GpEp=Y$Gp-Y$Ep, GpE=Y$Gp-Y$E)
  
  List_names_tot <- colnames(Z)
  
  Z_trunc <- Z 
  Z_trunc[Z_trunc==-2]=-1
  Z_trunc[Z_trunc==2]=1
  
  Y <- Y[keep_config,]
  Z <- Z[keep_config,]
  Z_trunc <- Z_trunc[keep_config,]
  num_conf <- 1:length(keep_config)
  row.names(Z_trunc) <- num_conf
  
  prob <- prob[keep_config]
  prob <- prob/sum(prob)
  
  # remove the duplicated due to the selection of variable 
  #+ 
  
  Y1 = Y
  rownames(Y1) = NULL
  #+ 
  
  remove_config = NULL
  
  if(!is.null(remove_var)){
    rm_ind = which(list_name_test == remove_var)
    Z_trunc_rm = Z_trunc[,-rm_ind]
    List_names_final = list_name_test[-rm_ind]
    rg_duplicate = which(duplicated2(Z_trunc_rm)) 
    keep_ind = duplicated3(Z_trunc_rm, Y1)
    remove_config = rg_duplicate[which(rg_duplicate %in% keep_ind == FALSE)]
  }
 
  Z_trunc_final <- Z_trunc[-remove_config,]
  
  prob_final <- prob[-remove_config] 
  prob_final <- prob_final/sum(prob_final)
  
  Z_trunc_final_code <- Z_trunc_final
  Z_trunc_final_code[Z_trunc_final_code==1]=3
  Z_trunc_final_code[Z_trunc_final_code==-1]=2
  Z_trunc_final_code[Z_trunc_final_code==0]=1
  # save the list of configurations
  head(Z_trunc_final_code)
  config_list <- 1:38
  config_list_final <- config_list[-remove_config] 
  saveRDS(Z_trunc_final_code , file = paste0(file_path_Results,"List_config.rds"))
}