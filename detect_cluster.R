#' Title
#'
#' @param noise_model 6 noise model to test significance 
#' @param Six_Series Data 
#' @param List_changes list of change-points including the begning and the end
#'
#' @return List of change-points after remove clusters
#' @return List of rejection day
#' @export 
#'
#' @examples

main_st = "0fal"
nearby_st = "0hoa"
seg_main = date_mean %>% filter(name == main_st)
seg_nb = date_mean %>% filter(name == nearby_st)
List_changes = list(main = c(seg_main$begin[1], seg_main$end),
                    nearby = c(seg_nb$begin[1], seg_nb$end)) 
Six_Series = read_data_new(path_data = path_data_NGL,
                           main_st = main_st,
                           nearby_st = nearby_st,
                           name_six_diff = c("GE", "GGp", "GEp", "EEp", "GpEp", "GpE")) 

remove_cluster <- function(noise_model, Six_Series, List_changes){
  
  Reject <- data.frame(Begin = as.Date(character()), 
                       End = as.Date(character()))
  Update_List_Changes <- as.Date(character(0))

  main_df = Six_Series %>% 
    select(Date, GE) %>% 
    drop_na()
  
  nearby_df = Six_Series %>% 
    select(Date, GpEp) %>% 
    drop_na()
  
  list_changes = List_changes$main
  length_seg <- numeric(length(list_changes) - 1)
  length_seg[1] <- sum(main_df$Date >= list_changes[1] &
                              main_df$Date <= list_changes[2])
  
  for (i in 2:(length(list_changes) - 1)) {
    length_seg[i] <- sum(main_df$Date > list_changes[i] &
                                main_df$Date <= list_changes[i + 1])
  }
  
  get_cluster <- function(X){
    output <- integer(length(X))
    counter <- 1 # Initialize a counter for numbering elements < 80
    
    for(i in 1:length(X)) {
      if(X[i] > 80) {
        output[i] <- 0
        counter <- 1 # Reset counter for the next cluster of elements < 80
      } else {
        output[i] <- counter
        counter <- counter + 1 # Increment counter for the next element in the current cluster
      }
    }
    
    return(output)
  }
  
  check_cluster = get_cluster(length_seg)
  
  if (all(check_cluster == 0)) {
    Reject = Reject
    Update_List_Changes = List_changes
  } else{
    ind_cluster = which(check_cluster == 1)
    if(1 %in% ind_cluster){
      # remove all changepoint in the clusters close to the begining 
      
    }
        
  }

}