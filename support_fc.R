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

test1 <- function(main_brp, nearby_brps){
  min1 = main_brp - 10
  max1 = main_brp + 10
  crenel_pos = NA
  crenel = nearby_brps[nearby_brps > min1 & nearby_brps < max1]
  if(length(crenel) > 0){
    d = crenel - main_brp
    crenel_pos = d[which.max(abs(d))]
  }
  
  return(crenel_pos)
}

test2 <- function(main_brp, nearby_brp){
  min1 = main_brp - 10
  max1 = main_brp + 10
  crenel_length = NA
  crenel = nearby_brp[nearby_brp > min1 & nearby_brp < max1]
  if(length(crenel) > 0){
    crenel_length = max(abs(crenel - main_brp))
  }
  return(crenel_length)  
}