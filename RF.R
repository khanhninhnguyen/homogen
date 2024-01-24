#' Train classifier RANDOM FOREST
#' INPUT:
#' * test results from FGLS 
#' * segmentation results of main and nearby station
#'  (date_mean file from segmentation.R) 
#' 
# Final predictive rule prog
library(caret)


# output: last predtive rule and prediction results 
# input: test result, name of series used for training: 
# G-E, G-G', G-E',E-E', G'-E', G'-E, remove them by name of remove variable 
# prob for each configuration, significance.level, name.version

path_code_att = path_code
path_restest <- paste0(path_results,"attribution/predictive_rule/")
file_path_Results=paste0(path_results,'attribution/predictive_rule/')
source(file = paste0(path_code, "support_RF.R"))

name.version ="FGLS_jump_tvalue.txt"
name.results <- paste0(path_restest, name.version) # name of test result file

# parameters 
# NbSim = 37865
significance.level = 0.05
B = 20
offset=0
GE=0
number.pop = 3

# make the truth table and probability ------------------------------------

G=c(rep(1,9), rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9))
E=c(rep(0,9),rep(-1,9),rep(0,9),rep(1,9),rep(-1,9),rep(1,9))
a=c(rep(0,3),rep(1,3),rep(-1,3))
Gp=rep(a,6)
Ep=rep(c(0,1,-1),18)
Y=data.frame(G=G,E=E,Gp=Gp,Ep=Ep)

Z=data.frame(GE=Y$G-Y$E,GGp=Y$G-Y$Gp,GEp=Y$G-Y$Ep,EEp=Y$E-Y$Ep,GpEp=Y$Gp-Y$Ep,GpE=Y$Gp-Y$E)
knitr::kable(head(Z))
List.names.tot <- colnames(Z)

Z.trunc <- Z 
Z.trunc[Z.trunc==-2]=-1
Z.trunc[Z.trunc==2]=1
knitr::kable(head(Z.trunc))

prob <- c(0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.010125,0.18225,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,0.010125,
          0.18225,0.010125,0.010125,0.010125,0.0005625,0.0005625,0.010125,0.0005625,0.0005625,
          0.010125,0.18225,0.010125,0.0005625,0.010125,0.0005625,0.0005625,0.010125,0.0005625,
          0.00225,0.00225,0.0405,0.000125,0.000125,0.00225,0.000125,0.000125,0.00225,0.00225,
          0.0405,0.00225,0.000125,0.00225,0.000125,0.000125,0.00225,0.000125)
keep.config <- c(1:3,6:15,17,19:24,26,28:30,33:40,43,46:49,52)
Y <- Y[keep.config,]
Z <- Z[keep.config,]
Z.trunc <- Z.trunc[keep.config,]
num.conf <- 1:length(keep.config)
row.names(Z.trunc) <- num.conf

prob <- prob[keep.config]
prob <- prob/sum(prob)

# remove the duplicated due to the selection of variable ------------------
# keep the one with break in GPS
#+ 
remove.var = "G-E"

Y1 = Y
rownames(Y1) = NULL
#+ 
rm.ind = which(list_name_test == remove.var)
Z.trunc.rm = Z.trunc[,-rm.ind]
List.names.final = list_name_test[-rm.ind]
sum(duplicated2(Z.trunc[,-rm.ind]))/2
rg.duplicate <- which(duplicated2(Z.trunc.rm)) 
keep.ind = duplicated3(Z.trunc.rm, Y1)

remove.config = rg.duplicate[which(rg.duplicate %in% keep.ind == FALSE)]
Z.trunc.final <- Z.trunc[-remove.config,]

prob.final <- prob[-remove.config] 
prob.final <- prob.final/sum(prob.final)

Z.trunc.final.code <- Z.trunc.final
Z.trunc.final.code[Z.trunc.final.code==1]=3
Z.trunc.final.code[Z.trunc.final.code==-1]=2
Z.trunc.final.code[Z.trunc.final.code==0]=1
# save the list of configurations
head(Z.trunc.final.code)
config.list <- 1:38
config.list.final <- config.list[-remove.config] 
saveRDS(Z.trunc.final.code , file = paste0(file_path_Results,"List_config.rds"))

#
List_main = read.table(file = paste0(path_restest,"list_selected_nmin200_10nearby.txt"), 
                       header = TRUE, stringsAsFactors = FALSE) 

Data.Res.Test0 <- read.table(name.results,
                             header = TRUE,
                             stringsAsFactors = FALSE) 

Data.Res.Test <- cbind(List_main[,c("main", "brp", "nearby")],
                       Data.Res.Test0[,7:12]) %>%
  mutate(brp = as.Date(List_main$brp, format="%Y-%m-%d")) 
colnames(Data.Res.Test)[4:9] <- paste0("t", List.names.tot) 


# Separate into learning and test dataset --------------------------------------------------------

NumDetected.per.station <- Data.Res.Test %>% group_by(main) %>% 
  summarise(ldetected.per.station=length(unique(brp))) %>% dplyr::select(ldetected.per.station)

NumNearby.per.station.detection <- Data.Res.Test %>% group_by(main,brp) %>% 
  count() %>% as.data.frame() %>% dplyr::select(n)

List.names.final = List.names.tot[-rm.ind]
for (i in  1:length(List.names.final)){
  eval(parse(text=paste0("Data.",List.names.final[i],"=Keep.Data(List.names.final[i])")))  
}

# equal prob data  --------------------------------------------------------
R = 100
Nbconfig <- nrow(Z.trunc.final)

error.test.4.methods <- matrix(NA,nrow=B,ncol=4)
NbSim <- R*Nbconfig #in order to have at least 5 samples for each configurations. 
# B <- 5
# NbSimLearn <- NbSim*0.8
# NbSimTest <- NbSim*0.2
set.seed(1)

error.test <- rep(NA, B)

for (b in 1:B){
  
  ######
  # Vfold 
  existing.pop.Learn <- 0
  existing.pop.Test <- 0
  
  while ((existing.pop.Learn!=15) & (existing.pop.Test !=15)){
    trainIndex<-createDataPartition(1:nrow(Data.GGp) , p=0.8, list = FALSE)
    
    Data.GGp.Learn <- Data.GGp %>% slice(trainIndex)
    Data.GGp.Test <- Data.GGp %>% slice(-trainIndex)
    
    Data.EEp.Learn <- Data.EEp %>% slice(trainIndex)
    Data.EEp.Test <- Data.EEp %>% slice(-trainIndex)
    
    Data.GEp.Learn <- Data.GEp %>% slice(trainIndex)
    Data.GEp.Test <- Data.GEp %>% slice(-trainIndex)
    
    Data.GpEp.Learn <- Data.GpEp %>% slice(trainIndex)
    Data.GpEp.Test <- Data.GpEp %>% slice(-trainIndex)
    
    Data.GpE.Learn <- Data.GpE %>% slice(trainIndex)
    Data.GpE.Test <- Data.GpE %>% slice(-trainIndex)
    
    
    existing.pop.Learn=sum(c(length(unique(Data.GGp.Learn$pop)),length(unique(Data.GEp.Learn$pop)),length(unique(Data.EEp.Learn$pop)),length(unique(Data.GpEp.Learn$pop)),length(unique(Data.GpE.Learn$pop)))) 
    
    existing.pop.Test=sum(c(length(unique(Data.GGp.Test$pop)),length(unique(Data.GEp.Test$pop)),length(unique(Data.EEp.Test$pop)),length(unique(Data.GpEp.Test$pop)),length(unique(Data.GpE.Test$pop))))  
  }
  
  #####
  # Construction of the pop of t.values for each series Learn and Test
  for (i in 1:length(List.names.final)){
    eval(parse(text=paste0("Pop.",List.names.final[i],".Learn=Pop.create(Data.",List.names.final[i],".Learn)")))
    eval(parse(text=paste0("Pop.",List.names.final[i],".Test=Pop.create(Data.",List.names.final[i],".Test)")))
  }
  
  
  #####
  # Bootstrap: ccontruction de DataLearn et DataTest respecting the proportion of configuration use Boot1, 
  # Boot() for the same number of configuration
  DataLearn <- c()
  DataTest <- c()
  DataLearn <- Boot("Learn",Z.trunc.final.code,NbSim) ####BOOT FC HAS PB
  DataTest <- Boot("Test",Z.trunc.final.code,NbSim)
  saveRDS(DataLearn, file = paste0(file_path_Results,"DataLearn_",b,significance.level, offset, GE, number.pop,".rds"))
  saveRDS(DataTest, file = paste0(file_path_Results,"DataTest_",b,significance.level, offset, GE, number.pop,".rds"))
  
  Res.pred1 <- PredRule_RF(DataLearn,DataTest,b,Nbconfig)
  saveRDS(Res.pred1, file = paste0(file_path_Results,"Res.pred_",b,significance.level, offset, GE, number.pop,".rds"))
  error.test[b] <- Res.pred1$err.tot 
  print(b)
}


