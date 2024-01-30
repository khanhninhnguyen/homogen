#' functions used in RF1
#' Construct the logical table 
#' @param prob probability of all configuration
#' @param keep_config index of configuration that want to be kept

# print all index where the value is duplicated
duplicated2 <- function(x){
  auxi <- do.call("paste",c(x,sep="_"))
  y <- table(auxi)
  y <- y > 1
  res <- y[match(auxi,names(y))]
  dimnames(res) <- NULL
  return(res)
}
# print index are duplicated but need to keep 1 for each pair,
# keep the one have G!=0
duplicated3 <- function(X, Y){
  auxi <- do.call("paste",c(X,sep="_"))
  ind.dup = which(duplicated(auxi))
  all.keep = sapply(c(1:length(ind.dup)), function(z){
    ind.dupli = which(auxi == auxi[ind.dup[z]])
    keep.ind = ind.dupli[which(Y[ind.dupli,1] !=0)]
    return(keep.ind)
  })
  
  return(all.keep)
}

construct_logical_table <- function(prob, keep_config, remove_var, list_name_test, path_save){
  G = rep(c(1, 0, -1, 0, 1, -1), each = 9)
  E = rep(c(0, -1, 0, 1, -1, 1), each = 9)
  Gp = rep(c(rep(0,3), rep(1,3), rep(-1,3)), times = 6)
  Ep = rep(c(0,1,-1), times = 18)
  
  Y <- data.frame(G, E, Gp, Ep)
  Z = data.frame(GE = Y$G - Y$E, GGp = Y$G - Y$Gp, GEp = Y$G - Y$Ep,
                 EEp = Y$E - Y$Ep, GpEp = Y$Gp - Y$Ep, GpE = Y$Gp - Y$E)
  
  List_names_tot <- colnames(Z)
  
  Z_trunc <- apply(Z, c(1,2), function(x) ifelse(x == -2 | x == 2, sign(x), x)) %>% 
    as.data.frame()
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
 
  Z_trunc_final <- Z_trunc[!seq_len(nrow(Z_trunc)) %in% remove_config, ]
  
  prob_final <- prob[!seq_along(prob) %in% remove_config]
  prob_final <- prob_final/sum(prob_final)
  
  Z_trunc_final_code <- Z_trunc_final
  Z_trunc_final_code[Z_trunc_final_code==1]=3
  Z_trunc_final_code[Z_trunc_final_code==-1]=2
  Z_trunc_final_code[Z_trunc_final_code==0]=1
  # save the list of configurations
  # head(Z_trunc_final_code)
  config_list <- 1:nrow(Z_trunc)
  config_list_final <- config_list[!seq_along(config_list) %in% remove_config]
  
  saveRDS(Z_trunc_final_code , file = paste0(path_save,"List_config.rds"))
  
  return(Z_trunc_final_code)
}

