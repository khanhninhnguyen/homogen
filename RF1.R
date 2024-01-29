#' Train classifier RANDOM FOREST
#' input:
#' * test results from FGLS 
#' output:
#' list of changepoints attributed 
#' Final predictive rule program 
#' 
#' 
####### Test result - remove when finish ##################
significance.level = 0.05
B = 20
offset=0
GE=0
number.pop = 3

library(caret)