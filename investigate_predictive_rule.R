# study all iterations giving 0 error 

all_pred = combine_rules(version = "ver4/", 
                         B,
                         file_path_Results,
                         Data_Res_Test,
                         significance_level,
                         offset,
                         GE, 
                         number_pop)
Final_table = read.table(file = paste0(file_path_Results, version = "ver4", "/FinalTable.txt"))
z_truth = Final_table$Z.truth

most_frequent_with_ties <- function(x) {
  uniq_vals <- unique(x)
  counts <- table(x)  # Get counts of each unique value
  max_counts <- max(counts)  # Find the highest frequency
  most_freq_vals <- names(counts)[counts == max_counts]  # Find all values with that frequency
  return(paste(most_freq_vals, collapse=", "))  # Return all such values, collapsed into a single string
}

max_freq <- function(x) {
  uniq_vals <- unique(x)
  counts <- table(x)  # Get counts of each unique value
  max_counts <- round(max(counts)*100/length(x), digits = 2)  # Find the highest frequency
  return(max_counts)  # Return all such values, collapsed into a single string
}

select_final_c <- function(x) {
  
  if (isTRUE(grepl(",",x))) {
    list_of_elements <- as.numeric(unlist((strsplit(x, ","))))
    prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
              0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
              0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
              0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
              0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
              0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
    y = list_of_elements[which.max(prob[list_of_elements])]
  } else {
    y = x
  }
  return(y)
} 


all_res = cbind(Data_Res_Test[,c(1:3)], all_pred) %>% 
  mutate(max_freq_config = apply(all_pred, 1, most_frequent_with_ties),
         freq = apply(all_pred, 1, max_freq),
         final_config = sapply(max_freq_config, select_final_c)) 
  
Final_table <- cbind(Final_table, all_res[,c(28:30)])
Final_table$pred.y <- NULL

write.table(Final_table, file = paste0(file_path_Results, 'ver4', "/FinalTable.txt"), sep = '\t', quote = FALSE)

# how many case have the same configuration? ----
check_same_value <- function(row) {
  as.numeric(length(unique(row)) == 1)
}

all_res$same_value <- apply(all_res[,4:15], 1, check_same_value)
table(all_res$same_value)

# can we decide when there are more than 1 predicted config? 
check_unique_value <- function(row) {
  c(unique(row))
}
check_unique_length <- function(row) {
  length(c(unique(row)))
}

check_freq_value <- function(row) {
  n = length(unique(row))
  l = rep(NA, n)
  for (i in c(1:n)) {
    l[i] <- round((length(which(row == unique(row)[i]))/length(row)), digits = 3)*100
  }
  return(l)
}

check_freq_diff <- function(row) {
  a = unlist(row)
  l = ifelse(length(unique(a))==1, 0,1)
  return(l)
}

max_freq_diff <- function(row) {
  a = unlist(row)
  l = max(a)
  return(l)
}


find_single_mode <- function(row) {
  freq <- table(row)
  return(names(freq)[which.max(freq)])
}

ind_diff = which(all_res$same_value !=1)
inves = all_res[ind_diff,]
# unique predicted config
inves$uniq <- apply(inves[,4:27], 1, check_unique_value)
inves$uniq_l <- apply(inves[,4:27], 1, check_unique_length)
# probability of predicted config
inves$freq <- apply(inves[,4:27], 1, check_freq_value)
# how many case have the same probability?
inves$diff<- sapply(inves$freq, FUN = check_freq_diff)
table(inves$diff)
a = inves %>% filter(diff==0)

inves$max_fre <- sapply(inves$freq, FUN = max_freq_diff)




# harmonize with Olivier results ------------------------------------------


