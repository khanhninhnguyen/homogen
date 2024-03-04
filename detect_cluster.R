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
  
  length_seg_main <- numeric(length(List_changes$main) - 1)
  length_seg_main[1] <- sum(main_df$Date >= List_changes$main[1] &
                              main_df$Date <= List_changes$main[2])
  
  for (i in 2:(length(List_changes$main) - 1)) {
    length_seg_main[i] <- sum(main_df$Date > List_changes$main[i] &
                                main_df$Date <= List_changes$main[i + 1])
  }
  
  length_seg_nearby <- numeric(length(List_changes$nearby) - 1)
  length_seg_nearby[1] <- sum(nearby_df$Date >= List_changes$nearby[1] &
                                nearby_df$Date <= List_changes$nearby[2])
  
  for (i in 2:(length(List_changes$nearby) - 1)) {
    length_seg_nearby[i] <- sum(nearby_df$Date > List_changes$nearby[i] & 
                                  nearby_df$Date <= List_changes$nearby[i + 1])
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
  
  check_cluster_main = get_cluster(length_seg_main)
  check_cluster_nearby = get_cluster(length_seg_nearby)
  check_cluster = check_cluster_main
  
  ind_cluster_main = which(check_cluster_main==1)
  ind_cluster_nearby = which(check_cluster_nearby==1)
  
  ind_cluster = which(check_cluster_main==1)
  if (all(check_cluster == 0)) {
    Reject = Reject
    Update_List_Changes = List_changes
  } else{
    
  }

}