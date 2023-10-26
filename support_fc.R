#' additional functions in change-points attribution.

extract_list_brp <- function(date_time_list){
  
  date_time_list = date_time_list %>% 
    group_by(name) %>%
    mutate(StationCount = n())
  
  filtered_list = filter(date_time_list, StationCount>1) 
  
  list_brp = filtered_list %>% 
    group_by(name) %>%
    mutate(Sequence = row_number()) %>%
    filter(Sequence!=1) %>%
    mutate(brp = as.Date(begin, format = "%Y-%m-%d"))
  
  List_brp = list_brp[, c("name", "brp")]
  
  return(List_brp)

}

read_data_new <- function(path_data, main_st, nearby_st, name_six_diff){
  # order: GPS-ERA, GPS-GPS', GPS-ERA', ERA-ERA', GPS'-ERA', GPS'-ERA
  
  name_pre = paste0("diwv_", main_st, "-", nearby_st, "_Zc")
  
  list_six_data = lapply(c(1:6), function(i){
    read.table(file = paste0(path_data, name_pre, i, ".txt"), 
               skip = 1, colClasses = c("character", "numeric"), 
               col.names = c("Date", name_six_diff[i]))
  })
  
  df_six_data = reduce(list_six_data, full_join, by = "Date")
  df_six_data$Date = convert_date(df_six_data$Date)
  
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

test2 <- function(main_brp, main_brps, df_data, nearby_brps, test_1){
  
  crenel_length_max = 62

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
    if(test_1 < 0){
      main_brp_min1 = main_brp + test_1 - 1
    }else if(test_1 > 0){
      main_brp_max1 = main_brp + test_1 + 1 
    }
  }
  
  # look for 2 closest brps: remove other changepoints in 62 days 
  
  other_dist = main_brps - main_brp
  noise = which(abs(other_dist) > 0 & abs(other_dist) < crenel_length_max)
  
  if(length(noise)>0){
    dist_noise = other_dist[noise]
    dist_noise = dist_noise[which.max(abs(dist_noise))]
  }else{
    dist_noise = 0
  }
  main_brp_min = main_brp - abs(dist_noise)
  main_brp_max = main_brp + abs(dist_noise)
  
  main_beg_new <- max(main_brps_m[main_brps_m < main_brp_min])
  main_end_new <- min(main_brps_m[main_brps_m > main_brp_max])
  
  list_nearby_brpb <- nearby_brps_m[nearby_brps_m < main_brp_min1]
  if(length(list_nearby_brpb) > 0){
    nearby_beg_new <- max(list_nearby_brpb)
  }else{
    nearby_beg_new <- 0
  }
  
  list_nearby_brpa <- nearby_brps_m[nearby_brps_m > main_brp_max1]
  if(length(list_nearby_brpa) > 0){
    nearby_end_new <- min(list_nearby_brpa)
  }else{
    nearby_end_new <- 0
  }
  
  out <- list(dist_noise = dist_noise, 
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
    remove_ind = which(df_data$Date < min(main_brp, main_brp + test_1) & 
                        df_data$Date > max(main_brp, main_brp + test_1))
    df_data[remove_ind, !(names(df_data) == "GPS_ERA")] <- NA
  }
  
  # check crenel in the main 
  
  if(dist_noise != 0){
    remove_ind = which(df_data$Date < min(main_brp, main_brp + dist_noise) & 
                         df_data$Date > max(main_brp, main_brp + dist_noise))
    df_data[remove_ind, !(names(df_data) == "GPS1_ERA1")] <- NA
  }
  
  df1_bef = df_data[which(df_data$Date > main_beg & df_data$Date < main_brp),]
  df1_aft = df_data[which(df_data$Date >= main_brp & df_data$Date < main_end),]
  
  df2_bef = df_data[which(df_data$Date > nearby_beg & df_data$Date < main_brp),]
  df2_aft = df_data[which(df_data$Date >= main_brp & df_data$Date < nearby_end),]
  
  df3_bef = df_data[which(df_data$Date > beg_df & df_data$Date < main_brp),]
  df3_aft = df_data[which(df_data$Date >= main_brp & df_data$Date < end_df),]
  
  if(all(is.na(df3_bef$GPS_GPS1))){
    nb3_bef = 0
  }else{
    nb3_bef = nb_consecutive(list.day = df3_bef$Date, x = df3_bef[["GPS1_ERA"]])
  }
  
  if(all(is.na(df3_aft$GPS_GPS1))){
    nb3_aft = 0
  }else{
    nb3_aft = nb_consecutive(list.day = df3_aft$Date, x = df3_aft[["GPS1_ERA"]])
  }
  
  if(all(is.na(df2_bef$GPS1_ERA1))){
    nb2_bef = 0
  }else{
    nb2_bef = nb_consecutive(list.day = df2_bef$Date, x = df2_bef[["GPS1_ERA1"]])
  }
  
  if(all(is.na(df2_aft$GPS1_ERA1))){
    nb2_aft = 0
  }else{
    nb2_aft = nb_consecutive(list.day = df2_aft$Date, x = df2_aft[["GPS1_ERA1"]])
  }
  
  nb1_bef = nb_consecutive(list.day = df1_bef$Date, x = df1_bef[["GPS_ERA"]])

  nb1_aft = nb_consecutive(list.day = df1_aft$Date, x = df1_aft[["GPS_ERA"]])

  period1_aft = as.numeric(max(df1_aft$Date) - min(df1_aft$Date))
  period2_aft = as.numeric(max(df2_aft$Date) - min(df2_aft$Date))
  period3_aft = as.numeric(max(df3_aft$Date) - min(df3_aft$Date))
  
  period1_bef = as.numeric(max(df1_bef$Date) - min(df1_bef$Date))
  period2_bef = as.numeric(max(df2_bef$Date) - min(df2_bef$Date))
  period3_bef = as.numeric(max(df3_bef$Date) - min(df3_bef$Date))
  
  r1_bef = nb1_bef/(period1_bef-1)
  r2_bef = nb2_bef/(period2_bef-1)
  r3_bef = nb3_bef/(period3_bef-1)
  
  r1_aft = nb1_aft/(period1_aft-1)
  r2_aft = nb2_aft/(period2_aft-1)
  r3_aft = nb3_aft/(period3_aft-1)
  
  return(list(
    n_main_bef = nb1_bef, n_nearby_bef = nb2_bef, n_joint_bef = nb3_bef,
    n_main_aft = nb1_aft, n_nearby_aft = nb2_aft, n_joint_aft = nb3_aft,
    r_main_bef = r1_bef, r_nearby_bef = r2_bef, r_joint_bef = r3_bef,
    r_main_aft = r1_aft, r_nearby_aft = r2_aft, r_joint_aft = r3_aft))
}