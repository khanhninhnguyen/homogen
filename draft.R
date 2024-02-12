version = "original_R100"
a = rep(NA,20)
for (b in c(1:20)) {
  FinalPred <- readRDS(paste0(file_path_Results,
                              version,
                              "/Res.pred_",
                              b,
                              significance_level, 
                              offset, 
                              GE, 
                              number_pop,
                              ".rds"))
  a[b] = FinalPred$err.tot
}
a
which.min(a)

# boxplot of the t value in the learning set 

learn_dat = a = readRDS(file = paste0(
  file_path_Results, version,
  "/DataLearn_", 
  b=1,
  significance_level, 
  offset, 
  GE, 
  number_pop,
  ".rds"))

boxplot(learn_dat[(80*20+1):(80*21),-6],
        main = paste0(version," configuration 22"),
        outline = FALSE)

boxplot(learn_dat[1:80,-6],
        main = paste0(version," configuration 1"),
        outline = FALSE)


boxplot(learn_dat[(320*20+1):(320*21),-6],
        main = paste0(version," configuration 22"),
        outline = FALSE)

boxplot(learn_dat[1:320,-6],
        main = paste0(version," configuration 1"),
        outline = FALSE)

dat = ori[which(ori$pred.y==22 & ver1$pred.y ==1),]
boxplot(dat[,5:9],
        main = paste0(version,"real"),
        outline = FALSE)

