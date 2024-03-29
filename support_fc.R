#' additional functions in change-points attribution.
get_clusters <- function(dates, threshold) {
  # dates must be sorted upward before
  if(length(dates) == 1){
    clusters_ind <- 0
  } else{
    clusters_ind <- rep(0, length(dates))
    sorted_dates <- sort(dates)
    clusters <- list()
    current_cluster <- c()
    
    for (i in 1:length(sorted_dates)) {
      if (length(current_cluster) == 0) {
        current_cluster <- c(sorted_dates[i])
      } else if (difftime(sorted_dates[i], tail(current_cluster, n = 1), units = "days") < threshold) {
        current_cluster <- c(current_cluster, sorted_dates[i])
      } else {
        clusters <- c(clusters, list(current_cluster))
        current_cluster <- c(sorted_dates[i])
      }
    }
    if (length(current_cluster) > 0) {
      clusters <- c(clusters, list(current_cluster))
    }
    list_clusters = clusters %>%
      keep(~ length(.) > 1)
    # rename clusters_ind 
    if(length(list_clusters) > 0){
      for (j in (1:length(list_clusters))) {
        clusterj = sort(list_clusters[[j]])
        for (k in c(1:length(clusterj))) {
          clusters_ind[which(dates == clusterj[k])] <- k
        }
      }
    }
  }

  return(clusters_ind)
}

extract_list_brp <- function(date_time_list){
  
  date_time_list = date_time_list %>% 
    group_by(name) %>%
    mutate(StationCount = n())
  
  filtered_list = filter(date_time_list, StationCount>1) 
  
  # List breaks from segmentation
  list_brp = filtered_list %>% 
    group_by(name) %>%
    mutate(Sequence = row_number()) %>%
    filter(Sequence!=max(StationCount)) %>%
    mutate(brp = as.Date(end, format = "%Y-%m-%d"))
  
  # remove the clusters of breaks by the first point
  
  List_brp = list_brp[, c("name", "brp")]
  
  return(List_brp)

}

read_data_new <- function(path_data, main_st, nearby_st, name_six_diff){
  # order: GPS-ERA, GPS-GPS', GPS-ERA', ERA-ERA', GPS'-ERA', GPS'-ERA
  
  name_pre = paste0("diwv_", main_st, "-", nearby_st, "_Zc")
  
  list_six_data = lapply(c(1:6), function(i){
    a = read.table(file = paste0(path_data, name_pre, i, ".txt"), 
               skip = 1, colClasses = c("character", "numeric"), 
               col.names = c("Date", name_six_diff[i])) %>%
      mutate(Date = as.Date(convert_date(Date), format = "%Y-%m-%d"))
  })
  
  df_six_data = reduce(list_six_data, full_join, by = "Date") %>%
    arrange(Date)
  
  return(df_six_data)
}

convert_date <- function(date_vec){
  date_out <- date_vec
  ind <- which(substr(date_vec, 1, 1) == "9")
  if (length(ind) == 0) {
    date_vec <- paste0("20", date_vec)
  } else {
    date_vec[ind] <- paste0("19", date_vec[ind])
    date_vec[-ind] <- paste0("20", date_vec[-ind])
  }
  
  date_out <- strptime(date_vec, "%Y%m%d.%H", tz = "GMT")
  date_out <- format(date_out, format = "%Y-%m-%d")
  class(date_out)
  return(date_out)
}

nb_consecutive <- function(list.day, x){
  a = list.day[which(is.na(x)== FALSE)]
  b = ts(a) 
  y = length(which(diff(b)==1))
  return(y)
}

test1 <- function(main_brp, nearby_brps){
  min1 = main_brp - 10
  max1 = main_brp + 10
  crenel_pos = NA
  crenel = nearby_brps[nearby_brps >= min1 & nearby_brps <= max1]
  if(length(crenel) > 0){
    d = crenel - main_brp
    crenel_pos = d[which.max(abs(d))]
  }
  
  return(crenel_pos)
}

test2 <- function(main_brp, main_brps, df_data, nearby_brps, test_1, cluster_ind){
  
  main_ind = which(!is.na(df_data$GPS_ERA))
  nearby_ind = which(!is.na(df_data$GPS1_ERA1))
  
  main_beg = df_data$Date[main_ind[1]]
  main_end = df_data$Date[main_ind[length(main_ind)]]
  
  nearby_beg = df_data$Date[nearby_ind[1]]
  nearby_end = df_data$Date[nearby_ind[length(nearby_ind)]]
  
  main_brps_m = c(main_brps,
                  main_beg, main_end)
  nearby_brps_m = c(nearby_brps, 
                    nearby_beg, nearby_end)
  
  # remove brpd which is coincident to main
  main_brp_min1 = main_brp
  main_brp_max1 = main_brp
  
  if(!is.na(test_1)){
    # if(test_1 < 0){
    #   main_brp_min1 = main_brp + test_1 - 1
    # }else if(test_1 > 0){
    #   main_brp_max1 = main_brp + test_1 + 1 
    # }
    main_brp_min1 = main_brp - abs(as.numeric(test_1)) - 1
    main_brp_max1 = main_brp + abs(as.numeric(test_1)) + 1 
  }
  
  main_brp_ind = which(main_brps == main_brp)
  noise = cluster_ind[main_brp_ind]
  if (noise == 1){
    noise_ind = which(main_brps == main_brp)
    diff_noise = diff(c(cluster_ind[-c(1:noise_ind)], 0))
    cluster_end_ind = which(diff_noise<0)[1]
    cluster_end = main_brps[(cluster_end_ind + main_brp_ind)]
    dist_noise = cluster_end - main_brp
  } else{
    dist_noise = 0
  }
  
  # limit cluster
  main_brp_min = main_brp - abs(dist_noise)
  main_brp_max = main_brp + abs(dist_noise)
  main_beg_new <- max(main_brps_m[main_brps_m < main_brp_min])
  if(main_beg_new != main_beg){
    main_beg_new <- main_beg_new + 1
  }
  main_end_new <- min(main_brps_m[main_brps_m > main_brp_max])

  list_nearby_brpb <- nearby_brps_m[nearby_brps_m < main_brp_min1]
  if(length(list_nearby_brpb) > 0){
    nearby_beg_new <- max(list_nearby_brpb)
    if(nearby_beg_new != nearby_beg)
    nearby_beg_new <- nearby_beg_new + 1
  }else{
    nearby_beg_new <- 0
  }
  
  list_nearby_brpa <- nearby_brps_m[nearby_brps_m > main_brp_max1]
  if(length(list_nearby_brpa) > 0){
    nearby_end_new <- min(list_nearby_brpa)
  }else{
    nearby_end_new <- 0
  }
  
  out <- list(dist_noise = dist_noise, noise = noise,
              main_beg_new = main_beg_new, main_end_new = main_end_new,
              nearby_beg_new = nearby_beg_new, nearby_end_new = nearby_end_new
  )
  return(out)
}

test3 <- function(main_brp, df_data, main_beg, main_end, nearby_beg, nearby_end,
                  test_1, dist_noise){

  beg_df = max(main_beg, nearby_beg)
  end_df = min(main_end, nearby_end)
  
  # remove points when brp in the main and nearby coincident
  if(!is.na(test_1)){
    remove_ind = which(df_data$Date < min(main_brp, main_brp - abs(as.numeric(test_1))) & 
                        df_data$Date > max(main_brp, main_brp + abs(as.numeric(test_1))))
    df_data[remove_ind, !(names(df_data) %in% c("Date", "GPS1_ERA1"))] <- NA
  }
  
  # check crenel in the main 
  
  if(dist_noise != 0){
    remove_ind = which(df_data$Date < min(main_brp, main_brp + dist_noise) & 
                         df_data$Date > max(main_brp, main_brp + dist_noise))
    df_data[remove_ind, !(names(df_data) %in% c("Date","GPS_ERA"))] <- NA
  }
  
  df1_bef = df_data[which(df_data$Date >= main_beg & df_data$Date <= main_brp),]
  df1_aft = df_data[which(df_data$Date > main_brp & df_data$Date <= main_end),]
  
  df2_bef = df_data[which(df_data$Date >= nearby_beg & df_data$Date <= main_brp),]
  df2_aft = df_data[which(df_data$Date > main_brp & df_data$Date <= nearby_end),]
  
  df3_bef = df_data[which(df_data$Date >= beg_df & df_data$Date <= main_brp),]
  df3_aft = df_data[which(df_data$Date > main_brp & df_data$Date <= end_df),]
  
  
  nb1_bef = length(which(!is.na(df1_bef[["GPS_ERA"]])))
  nb1_aft = length(which(!is.na(df1_aft[["GPS_ERA"]])))
  nb2_bef = length(which(!is.na(df2_bef[["GPS1_ERA1"]])))
  nb2_aft = length(which(!is.na(df2_aft[["GPS1_ERA1"]])))
  nb3_bef = length(which(!is.na(df3_bef[["GPS1_ERA"]])))
  nb3_aft = length(which(!is.na(df3_aft[["GPS1_ERA"]])))
  
  # if(all(is.na(df3_bef$GPS_GPS1))){
  #   nb3_bef = 0
  # }else{
  #   nb3_bef = nb_consecutive(list.day = df3_bef$Date, x = df3_bef[["GPS1_ERA"]])
  # }
  # 
  # if(all(is.na(df3_aft$GPS_GPS1))){
  #   nb3_aft = 0
  # }else{
  #   nb3_aft = nb_consecutive(list.day = df3_aft$Date, x = df3_aft[["GPS1_ERA"]])
  # }
  # 
  # if(all(is.na(df2_bef$GPS1_ERA1))){
  #   nb2_bef = 0
  # }else{
  #   nb2_bef = nb_consecutive(list.day = df2_bef$Date, x = df2_bef[["GPS1_ERA1"]])
  # }
  # 
  # if(all(is.na(df2_aft$GPS1_ERA1))){
  #   nb2_aft = 0
  # }else{
  #   nb2_aft = nb_consecutive(list.day = df2_aft$Date, x = df2_aft[["GPS1_ERA1"]])
  # }
  # 
  # nb1_bef = nb_consecutive(list.day = df1_bef$Date, x = df1_bef[["GPS_ERA"]])
  # nb1_aft = nb_consecutive(list.day = df1_aft$Date, x = df1_aft[["GPS_ERA"]])

  period1_aft = as.numeric(max(df1_aft$Date) - min(df1_aft$Date))
  period2_aft = as.numeric(max(df2_aft$Date) - min(df2_aft$Date))
  period3_aft = as.numeric(max(df3_aft$Date) - min(df3_aft$Date))

  period1_bef = as.numeric(max(df1_bef$Date) - min(df1_bef$Date))
  period2_bef = as.numeric(max(df2_bef$Date) - min(df2_bef$Date))
  period3_bef = as.numeric(max(df3_bef$Date) - min(df3_bef$Date))
  
  r1_bef = nb1_bef/(period1_bef)
  r2_bef = ifelse(period2_bef <= 1, 0, nb2_bef/(period2_bef))
  r3_bef = ifelse(period2_bef <= 1, 0, nb3_bef/(period3_bef))
  
  r1_aft = nb1_aft/(period1_aft-1)
  r2_aft = ifelse(period2_aft <= 1, 0, nb2_aft/(period2_aft))
  r3_aft = ifelse(period3_aft <= 1, 0, nb3_aft/(period3_aft))
  
  return(list(
    n_main_bef = nb1_bef, n_nearby_bef = nb2_bef, n_joint_bef = nb3_bef,
    n_main_aft = nb1_aft, n_nearby_aft = nb2_aft, n_joint_aft = nb3_aft,
    r_main_bef = r1_bef, r_nearby_bef = r2_bef, r_joint_bef = r3_bef,
    r_main_aft = r1_aft, r_nearby_aft = r2_aft, r_joint_aft = r3_aft))
}

extract_info_nearby <- function(path_data, list_brp, path_results){
  
  list_all_pairs = unique(substr(list.files(path_data),6,14))
  list_main = substr(list_all_pairs,1,4)
  list_nearby = substr(list_all_pairs,6,9)
  
  list_main_unique = unique(list_main)
  list_main_inhomo = list_main_unique[which(list_main_unique %in% list_brp$name)]
  
  out = data.frame(main = character(0),
                   brp = as.Date(character(0)),     
                   nearby = character(0),
                   coincident = numeric(0),
                   dist_noise = numeric(0),
                   noise = numeric(0),
                   main_beg_new = as.Date(character(0)),      
                   main_end_new = as.Date(character(0)),        
                   nearby_beg_new = as.Date(character(0)),     
                   nearby_end_new = as.Date(character(0)),      
                   n_main_bef = numeric(0),
                   n_nearby_bef = numeric(0),
                   n_joint_bef = numeric(0),
                   n_main_aft = numeric(0),
                   n_nearby_aft = numeric(0), 
                   n_joint_aft = numeric(0),
                   r_main_bef = numeric(0),
                   r_nearby_bef = numeric(0),
                   r_joint_bef = numeric(0),  
                   r_main_aft = numeric(0),
                   r_nearby_aft = numeric(0),   
                   r_joint_aft = numeric(0))
  
  for (i in c(1:length(list_main_inhomo))) {
    print(i)
    main_st = list_main_inhomo[i]
    main_brps = list_brp$brp[which(list_brp$name == main_st)]
    nearby_sts <- list_nearby[which(list_main == main_st)]
    cluster_indi = list_brp$cluster_index[which(list_brp$name == main_st)]
    
    for (j in c(1:length(main_brps))) {
      main_brp = main_brps[j]
      
      for (k in c(1:length(nearby_sts))) {
        nearby_st = nearby_sts[k]
        nearby_brps = list_brp$brp[which(list_brp$name == nearby_st)]
        # merge data 
        df_data = read_data_new(path_data = path_data,
                                main_st = main_st, 
                                nearby_st = nearby_st,
                                name_six_diff = name_six_diff)
        df_data$Date <- as.Date(df_data$Date, format = "%Y-%m-%d")
        # Test 1 
        
        test_1 <- test1(main_brp = main_brp, nearby_brps = nearby_brps)
        # Test 2 
        
        test_2 <- test2(main_brp = main_brp,
                        main_brps = main_brps, 
                        df_data = df_data,
                        nearby_brps = nearby_brps, 
                        test_1 = test_1,
                        cluster_ind = cluster_indi)
        if(test_2$nearby_beg_new !=0 & test_2$nearby_end_new != 0){
          test_3 <- test3(main_brp = main_brp, df_data = df_data,
                          test_1 = test_1, dist_noise = test_2$dist_noise,
                          main_beg = test_2$main_beg_new, 
                          main_end = test_2$main_end_new,
                          nearby_beg = test_2$nearby_beg_new,
                          nearby_end = test_2$nearby_end_new)
          
          out_ind = data.frame(main = main_st, brp = main_brp, nearby = nearby_st,
                               coincident = test_1)
          out_ind_all = out_ind %>% 
            bind_cols(as.data.frame(test_2)) %>% 
            bind_cols(as.data.frame(test_3))
          out <- out %>% 
            rbind(out_ind_all)
        }
      }
    }
  }
  out[,c(17:22)] <- out[,c(17:22)] %>%
    mutate(across(where(is.numeric), ~formatC(.x, format = "f", digits = 5)))
  
  write.table(out, file = paste0(path_results, "pre_info_test.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  
  return(out)
}

check_selected <- function(list_brp, infor_all, nbcsv_min, distance, rate_consecutive){
  filtered <- infor_all %>%
    filter(noise < 2, 
           n_main_bef > nbcsv_min, 
           n_nearby_bef > nbcsv_min, 
           n_joint_bef > nbcsv_min, 
           n_main_aft > nbcsv_min, 
           n_nearby_aft > nbcsv_min, 
           n_joint_aft > nbcsv_min,
           dd < distance,
           r_main_bef > rate_consecutive, 
           r_nearby_bef > rate_consecutive, 
           r_joint_bef > rate_consecutive, 
           r_main_aft > rate_consecutive, 
           r_nearby_aft > rate_consecutive, 
           r_joint_aft > rate_consecutive)  
  # check the number of change-points can be tested per stations 
  total_brp <- list_brp %>%
    group_by(name) %>%
    summarize(Count = n())
  
  filtered_brp <- unique( filtered[,c(1:2)]) %>%
    group_by(main) %>%
    summarize(Count = n())
  
  filtered_nb <- unique( filtered[,c(1:3)]) %>%
    group_by(main,brp) %>%
    summarize(Count = n())
  
  joint_df = left_join(filtered_brp, total_brp, by = join_by(main==name) ) %>%
    mutate(rate = Count.x/Count.y)
  
  out <- list(nb_main = length(unique(filtered_brp$main)),
              nb_nb = length(unique(filtered$nearby)),
              total_tested_brp = sum(joint_df$Count.x),
              total_brp = sum(joint_df$Count.y),
              nb_triplet = nrow(filtered), 
              full_main = nrow(joint_df[which(joint_df$rate==1),]),
              rate = summary(filtered_nb$Count))
  
  return(out)
}

list_longest_segment <- function(path_data, date_mean, list_brp, path_results){
  list_all_pairs = unique(substr(list.files(path_data),6,14))
  list_main = substr(list_all_pairs,1,4)
  list_nearby = substr(list_all_pairs,6,9)
  
  test_pairs <- unique(if_else(list_main<list_nearby, 
                               paste(list_main, list_nearby, sep = "-"), 
                               paste(list_nearby, list_main, sep = "-")))
  
  list_points = date_mean %>%
    group_by(name) %>%
    summarise(
      FirstBegin = first(begin),
      AllEnd = list(end)
    ) %>%
    ungroup() %>%
    split(.$name) %>%
    lapply(function(x) c(FirstBegin = x$FirstBegin, AllEnd = unlist(x$AllEnd)))
  
  df_infor <- as.data.frame(replicate(9, rep(NA, length(test_pairs)), 
                                      simplify = FALSE)) %>%
    mutate(name = test_pairs) %>%
    relocate(name, .before = 1) %>%
    setNames(c("name", paste0(rep(c("length_", "beg_", "end_"), each = 3),
                              c("main", "nearby", "joint")))) %>%
    mutate(across(5:10, \(x) as.Date(x, origin = "1970-01-01", format = "%Y-%m-%d"))) 
  
  for (i in c(1:length(test_pairs))) {
    station1 = substr(test_pairs[i], 1, 4)
    station2 = substr(test_pairs[i], 6, 9)
    
    list_brp1 = sort(as.Date(list_points[[which(names(list_points) == station1)]],
                             origin = "1970-01-01",
                             format = "%Y-%m-%d"))
    list_brp2 = sort(as.Date(list_points[[which(names(list_points) == station2)]],
                             origin = "1970-01-01",
                             format = "%Y-%m-%d"))
    
    beg = max(list_brp1[1], list_brp2[1])
    end = min(list_brp1[length(list_brp1)],list_brp2[length(list_brp2)])
    list_brp_all = sort(c(list_brp1, list_brp2))
    list_brp3 = list_brp_all[which(list_brp_all >= beg & list_brp_all <= end)]
    
    d1 = diff(list_brp1, lag = 1)
    d2 = diff(list_brp2, lag = 1)
    d3 = diff(list_brp3, lag = 1)
    
    df_infor[i,"beg_main"] = list_brp1[which.max(d1)]
    df_infor[i,"end_main"] = list_brp1[which.max(d1)+1]
    df_infor[i,"beg_nearby"] = list_brp2[which.max(d2)]
    df_infor[i,"end_nearby"] = list_brp2[which.max(d2)+1]
    df_infor[i,"beg_joint"] = list_brp3[which.max(d3)]
    df_infor[i,"end_joint"] = list_brp3[which.max(d3)+1]
    
    df_infor[i,paste0(rep("length_", 3), c("main", "nearby", "joint"))] =
      c(max(d1), max(d2), max(d3))
    
    # priotize the series has main station in list of main 
    two_stations = c(station1, station2) 
    ind_main = which(two_stations %in% list_main)
    if(2 %in% ind_main){
      add_row = df_infor[i,] %>% 
        mutate(name = paste(station2, station1, sep = "-"))
      add_row[,c("length_main", "length_nearby")] <- add_row[,c("length_nearby", "length_main")] 
      add_row[,c("beg_main", "beg_nearby")] <- add_row[,c("beg_nearby", "beg_main")] 
      add_row[,c("end_main", "end_nearby")] <- add_row[,c("end_nearby", "end_main")] 
      if(length(ind_main) == 2){
        df_infor = df_infor %>% bind_rows(add_row)
      }else if(ind_main == 2){
        df_infor[i,] <- add_row
      }
    }
  }
  
  
  df_infor <- df_infor %>%
    mutate(main = substr(name, 1, 4),
           nearby = substr(name, 6, 9)) %>%
    filter(name %in% list_all_pairs) %>%
    select(main, nearby, everything(), -name)
  
  write.table(df_infor, file = paste0(path_results, "List_longest_segment.txt"), 
              sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(df_infor)
}

# Define a custom function, respect the name of columns 
select_rows_based_on_conditions <- function(df) {
  df_filtered <- df %>% filter(dd < 50)
  
  if (nrow(df_filtered) > 10) {
    df_sorted <- df_filtered %>% 
      arrange(desc(n_joint_min)) %>%
      slice_max(order_by = n_joint_min, n = 10)
  } else if (nrow(df_filtered) < 10) {
    df_sorted <- df %>% 
      arrange(dd) %>% 
      slice_head(n = 10)
  } else {
    df_sorted <- df_filtered
  }
  
  return(df_sorted)
}

extract_FGLS_result <- function(list_selected_cases, path_FGLS){
  # standard output contains jumps and t-values
  # additional results: other estimates, ARMA coefficient, MW var 
  list_selected_cases = selected_cases
  nr <- nrow(list_selected_cases)
  
  main_res <- data.frame(matrix(NA, ncol = 12, nrow = nr))
  add_res_estimates <- data.frame(matrix(NA, ncol = 54, nrow = nr))
  add_res_arma <- data.frame(matrix(NA, ncol = 12, nrow = nr))
  add_res_MWvar <- data.frame(matrix(NA, ncol = 6, nrow = nr))
  
  for (i in c(1:nr)) {
    name_i = paste0(list_selected_cases$main[i], 
                    list_selected_cases$brp[i], 
                    list_selected_cases$nearby[i])
    test_i = get(load(paste0(path_FGLS, 
                             name_i,
                             "fgls.RData")))
    jump_est = sapply(c(1:6), function(x) test_i[[name_six_diff[x]]]$t.table$Estimate[9]) %>% 
      map_dbl(~ if(is.null(.x)) NA_real_ else .x)
    t_values = sapply(c(1:6), function(x) test_i[[name_six_diff[x]]]$t.table$`t value`[9]) %>% 
      map_dbl(~ if(is.null(.x)) NA_real_ else .x) 

    other_est = sapply(c(1:6), function(x) test_i[[name_six_diff[x]]]$t.table$Estimate[-9])
    if (is.null(other_est[[1]])) {
      other_est[[1]] <- rep(NA, 9)
    }
    
    arma_coef <- sapply(c(1:6), function(x) test_i[[name_six_diff[x]]]$coef.arma)
    if (is.null(arma_coef[[1]])) {
      arma_coef[[1]] <- rep(NA, 2)
    }
    
    MWvar <- sapply(1:6, function(x) {
      if (is.null(test_i[[name_six_diff[x]]])) {
        return(NA)  
      } else {
        return(mean(test_i[[name_six_diff[x]]]$var, na.rm = TRUE))  
      }
    })

    
    main_res[i,] = c(jump_est, t_values)
    add_res_estimates[i,] = unlist(other_est)
    add_res_arma[i,] = unlist(arma_coef)   
    add_res_MWvar[i,] = MWvar
    
  }
                       
  main_res %>% 
    mutate(across(everything(), ~sprintf("%.4f", .))) %>% 
    setNames(c(paste0("Jump_", name_six_diff), paste0("Tvalue_", name_six_diff))) %>%
    write.table(file = paste0(path_results, "FGLS_jump_tvalue.txt"), 
                sep = "\t", 
                row.names = FALSE, 
                col.names = TRUE, 
                quote = FALSE)
         
  add_res_estimates %>% 
    mutate(across(everything(), ~sprintf("%.4f", .))) %>% 
    setNames(as.vector(outer(
      rownames(test_i$GPS_GPS1$coefficients)[-9], 
      name_six_diff, paste, sep = "-"))) %>% 
    write.table(file = paste0(path_results, "FGLS_other_estimates.txt"), 
              sep="\t", 
              col.names = TRUE,
              row.names = FALSE, 
              quote = FALSE)
  
  add_res_arma %>% 
    mutate(across(everything(), ~sprintf("%.4f", .))) %>% 
    setNames( paste(rep(c("Phi", "Theta"),6), 
                    rep(name_six_diff,each =2), sep = "-")) %>% 
  write.table(file = paste0(path_results, "FGLS_arma_coef.txt"), 
              sep="\t", 
              col.names = TRUE,
              row.names = FALSE, 
              quote = FALSE)
  
  add_res_MWvar  %>% 
    mutate(across(everything(), ~sprintf("%.4f", .))) %>% 
    setNames(name_six_diff) %>% 
  write.table(file = paste0(path_results, "FGLS_MWvar.txt"), 
              sep="\t", 
              col.names = TRUE,
              row.names = FALSE, 
              quote = FALSE)
  
}


#' plot time series 

plot_test_res <- function(main_st, brp, nearby_st, 
                          name_nearby_full,
                          name_six_diff,
                          distance,
                          path_data_NGL, 
                          path_FGLS){
  
  test_case = get(load(paste0(path_FGLS, 
                              paste0(main_st, brp, nearby_st),
                              "fgls.RData")))
  test_main = get(load(paste0(path_FGLS, 
                              paste0(main_st, brp, name_nearby_full),
                              "fgls.RData")))
  
  df_data = read_data_new(path_data = path_data_NGL,
                          main_st = main_st, 
                          nearby_st = nearby_st,
                          name_six_diff = name_six_diff)
  
  six_vars = names(df_data)[!names(df_data) %in% c("Date")]
  four_vars <- six_vars[!six_vars %in% c("GPS_ERA", "GPS1_ERA1")]
  
  mean_vec = rep(NA, 12)
  mean_vec[c(1:2)] <- test_main$GPS_ERA$coefficients[c("left", "right"), ]
  mean_vec[3:12] <- unlist(lapply( six_vars[!six_vars %in% c("GPS_ERA")], function(x) test_case[[x]]$coefficients[c("left", "right"),]))
  
  main_beg = test_main$GPS_ERA$design.matrix$date[1]
  main_end = test_main$GPS_ERA$design.matrix$date[length(test_main$GPS_ERA$design.matrix$date)]
  nearby_beg = test_case$GPS1_ERA1$design.matrix$date[1]
  nearby_end = test_case$GPS1_ERA1$design.matrix$date[length(test_case$GPS1_ERA1$design.matrix$date)]
  joint_beg = test_case$GPS_GPS1$design.matrix$date[1] 
  joint_end = test_case$GPS_GPS1$design.matrix$date[length(test_case$GPS_GPS1$design.matrix$date)]
  
  df_data = df_data %>% 
    mutate(GPS_ERA = if_else(Date < main_beg | Date > main_end, NA_real_, GPS_ERA)) %>% 
    mutate(GPS1_ERA1 = if_else(Date < nearby_beg | Date > nearby_end, NA_real_, GPS1_ERA1)) %>% 
    mutate(across(all_of(four_vars),
                  ~if_else(Date < joint_beg | Date > joint_end, NA_real_, .)))
  
  
  df_data <- df_data[rowSums(is.na(df_data[, six_vars])) < length(six_vars), ] %>%
    complete(Date = seq(min(Date), max(Date), by = "day"))
  
  for (i in c(1:6)) {
    ind_brp = which(df_data$Date == brp)
    df_data[1:ind_brp, paste0("mean.", six_vars[i])] <- mean_vec[2*i-1]
    df_data[ind_brp:(nrow(df_data)), paste0("mean.", six_vars[i])] <- (mean_vec[2*i-1] + mean_vec[2*i])
    ind_na = which(is.na(df_data[six_vars[i]]))
    df_data[ind_na, paste0("mean.", six_vars[i])] <- NA
  }
  
  df_long <- pivot_longer(df_data,
                          cols = -Date,
                          names_to = "Variable", 
                          values_to = "Value") %>% 
    mutate(Type = ifelse(startsWith(Variable, "mean"), "Mean",  "IWVdiff")) %>% 
    mutate(Variable = gsub("mean\\.", "", Variable))

  # Custom labeller function
  my_labeller <- function(test_name) {
    if(test_name == "GPS_ERA"){
      Jump = test_main$GPS_ERA$coefficients[9]
      Tvalue = test_main$GPS_ERA$t.table$`t value`[9]
    }else{
      Jump = test_case[[test_name]]$coefficients[9]
      Tvalue = test_case[[test_name]]$t.table$`t value`[9]
      }
    return(paste0(test_name, 
                 " : Jump = ", sprintf("%.3f", Jump),
                 ", T-value = ", sprintf("%.3f", Tvalue),
                 ", Std. Err = ", sprintf("%.3f", Jump/Tvalue)))
  }
  
  df_long$FacetVar = sapply(df_long$Variable, my_labeller)

  FacetVar = as.character(unique(df_long$FacetVar))
  measurement_names <- sapply(strsplit(FacetVar, " :"), `[`, 1)
  order_index <- match(c("GPS_ERA", "ERA_ERA1", "GPS_ERA1", "GPS1_ERA1", "GPS_GPS1", "GPS1_ERA"), 
                       measurement_names)
  ordered_measurements <- FacetVar[order_index]
  
  df_long$FacetVar <- factor(df_long$FacetVar, 
                             levels = ordered_measurements)

  p <- ggplot(df_long, aes(x = Date, y = Value, color = Type)) + 
    theme_minimal() +
    geom_line(lwd = 0.2) +
    facet_wrap(~FacetVar, scales = "free_y", ncol = 2) +  # Using custom labeller
    geom_vline(xintercept = brp, lty = 3, lwd = 0.2) + 
    ylim(min(df_long$Value, na.rm = TRUE), max(df_long$Value, na.rm = TRUE)) + 
    scale_color_manual(values = c("IWVdiff" = "lightblue", "Mean" = "red")) +
    ylab(expression(paste("IWV Difference", (kg.m^-2)))) + 
    xlab("") + 
    theme(axis.text.x = element_text(size = 4.5), 
          axis.text.y = element_text(size = 5),
          legend.text = element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.2, "cm"), 
          plot.tag = element_text(size = 6) ,
          plot.subtitle = element_text(size = 6),
          legend.title = element_blank(), 
          legend.box.spacing = unit(0, "pt"),
          strip.text = element_text(size = 4), 
          panel.grid.major = element_line(size = 0.1, color = "grey80"),  # Default might be around 0.5
          panel.grid.minor = element_line(size = 0.05, color = "grey85"), # Default might be around 0.25
          plot.margin = rep(unit(0,"null"),4))
    
  jpeg(paste0(path_results,"figure/", main_st, brp, nearby_st, "dist", distance ,".jpeg"),
       width = 3000, height = 1800, res = 600) # change name
  
  print(p)
  # print(ml1)
  dev.off() 
  
}


plot_full_series <- function(main_st, nearby_st, 
                             path_data_NGL, 
                             date_mean,
                             name_six_diff,
                             distance){
  
  main_segments = date_mean %>% filter(name == main_st)
  nearby_segments = date_mean %>% filter(name == nearby_st)
  
  df_data = read_data_new(path_data = path_data_NGL,
                          main_st = main_st, 
                          nearby_st = nearby_st,
                          name_six_diff = name_six_diff) %>%
    tidyr::complete(Date = seq(min(Date), max(Date), by = "day")) %>%
    mutate(Mean_GPS_ERA = main_segments$mean[1],
           Mean_GPS1_ERA1 = nearby_segments$mean[1])
  
  if(nrow(main_segments) >1){
    for (j in 2:nrow(main_segments)) {
      ind_mean = which(df_data$Date <= main_segments$end[j] &
                         df_data$Date >= main_segments$begin[j])
      df_data$Mean_GPS_ERA[ind_mean] <- main_segments$mean[j] 
      df_data$Mean_GPS_ERA[which(is.na(df_data$GPS_ERA))] <- NA
    }
  }
  
  if(nrow(nearby_segments)>1){
    for (k in 2:nrow(main_segments)) {
      ind_mean = which(df_data$Date <= nearby_segments$end[k] &
                         df_data$Date >= nearby_segments$begin[k])
      df_data$Mean_GPS1_ERA1[ind_mean] <- nearby_segments$mean[k] 
      df_data$Mean_GPS1_ERA1[which(is.na(df_data$GPS1_ERA1))] <- NA
    }
  }
  
  plot_list <- list()
  
  for (i in c(1:6)) {
    # Access the design.matrix dataframe
    df <- df_data[, c("Date", name_six_diff[i])]
    colnames(df)[2] <- "Signal"
    
    ylab = name_six_diff[i]
    
    
    ymax = max(df_data[,-1], na.rm = TRUE)
    ymin = min(df_data[,-1], na.rm = TRUE)
  
    if (name_six_diff[i] %in% c("GPS_ERA", "GPS1_ERA1")){
      df <- df %>% mutate(Mean = df_data[[paste0("Mean_", name_six_diff[i])]])
      df_long <- pivot_longer(df, cols = -Date, names_to = "Variable", values_to = "Value")
      df_long$Variable <- factor(df_long$Variable, levels = c("Signal", "Mean"))
      if (name_six_diff[i] == "GPS_ERA"){
        brp = main_segments$end
      }else{
        brp = nearby_segments$end
      }      
      brp = as.Date(brp, format="%Y-%m-%d")
      brp = brp[-length(brp)]
      
      p <- ggplot(df_long, aes(x = Date, y = Value, col = Variable)) + 
        theme_minimal() + 
        geom_line(lwd = 0.2) + # or geom_point() depending on your data
        scale_color_manual(values = c("Signal" = "lightblue", "Mean" = "red")) +
        geom_vline(xintercept = brp, lty = 3, lwd = 0.2)  
    }else{
      p <- ggplot(df, aes(x = Date, y = Signal)) + 
        theme_minimal() + 
        geom_line(color = "lightblue", lwd = 0.2)
    }
    
    p <- p + 
      labs(subtitle = name_six_diff[i]) + 
      ylab(expression(paste("IWV Difference ", (kg.m^-2)))) + 
      xlab("") + 
      ylim(ymin, ymax) + 
      theme(
        axis.text.x = element_text(size = 4.5), 
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size=4),
        axis.title = element_text(size = 5), 
        legend.position = "none",
        plot.tag = element_text(size = 6) ,
        plot.subtitle = element_text(size = 6),
        legend.title = element_blank(), 
        legend.box.spacing = unit(0, "pt"),
        strip.text = element_text(size = 4), 
        panel.grid.major = element_line(size = 0.1, color = "grey80"),  # Default might be around 0.5
        panel.grid.minor = element_line(size = 0.05, color = "grey85"), # Default might be around 0.25
        plot.margin = rep(unit(0,"null"),4))
    # Add the plot to the list
    plot_list[[i]] <- p
  }
  jpeg(paste0(path_results,"figure/", main_st, nearby_st , "dist", distance, ".jpeg"),
       width = 3000, height = 1800, res = 600) # change name
  
  grid.arrange(grobs = plot_list, nrow = 3, ncol = 2)
  
  dev.off() 
}



