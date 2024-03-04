#' Title
#'
#' @param noise_model 
#' @param Six_Series 
#' @param List_changes 
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
  main_df = Six_Series %>% 
    select(Date, GE) %>% 
    drop_na()
  
  nearby_df = Six_Series %>% 
    select(Date, GpEp) %>% 
    drop_na()
  
  
  
}