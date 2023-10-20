#' additional functions in change-points attribution.

extract_list_brp <- function(date_time_list){
  
  date_time_list = date_time_list %>% 
    group_by(name) %>%
    mutate(StationCount = n())
  
  filtered_list = filter(date_time_list, StationCount>1) 
  
  list_brp = filtered_list %>% 
    group_by(name) %>%
    mutate(Sequence = row_number()) %>%
    filter(Sequence!=1) 
    
  List_brp = as.Date(list_brp$begin, format = "%Y-%m-%d") 
  
  return(List_brp)

}