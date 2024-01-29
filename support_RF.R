# all functions used for the predictive rule 
# print all index where the value is duplicated
duplicated2 <- function(x){
  auxi <- do.call("paste",c(x,sep="_"))
  y <- table(auxi)
  y <- y>1
  res <- y[match(auxi,names(y))]
  dimnames(res) <- NULL
  return(res)
}
# print index are duplicated but need to keep 1 for each pair,
# keep the one have G!=0
duplicated3 <- function(X, Y){
  auxi <- do.call("paste",c(X,sep="_"))
  ind.dup = which(duplicated(auxi))
  all.keep = sapply(c(1:length(ind.dup)), function(z){
    ind.dupli = which(auxi == auxi[ind.dup[z]])
    keep.ind = ind.dupli[which(Y[ind.dupli,1] !=0)]
    return(keep.ind)
  })
  
  return(all.keep)
}
# Separate data for each series 
Keep.Data <- function(name.series){
  Thresh <- significance.level
  names.t.series <- paste0("t",name.series)
  colnames.info <- c("main","brp","nearby")
  Data.Res <- Data.Res.Test[colnames(Data.Res.Test) %in% c(colnames.info,names.t.series)] 
  Data.Res$t <- Data.Res[,which(colnames(Data.Res) %in% names.t.series)]
  Data.Res <- Data.Res[,-which(colnames(Data.Res) %in% names.t.series)]
  Data.Res <- Data.Res %>% dplyr::mutate(p=2*pnorm(-abs(t))) %>%
    mutate(signif=ifelse(p<Thresh,1,0),
           t.sign=sign(t),
           pop=(signif==0)*1+((signif==1) & (t.sign==-1))*2+ ((signif==1) & (t.sign==1))*3,
           pop.code=(signif==0)*0+((signif==1) & (t.sign==-1))*-1+ ((signif==1) & (t.sign==1))*1)
  return(Data.Res)
}
# Creat population
Pop.create <- function(dataset){
  Pop1 <- dataset %>% dplyr::filter(pop==1) %>% dplyr::select(t)  %>% unlist()
  print(Pop1)
  Pop2 <- dataset %>% dplyr::filter(pop==2) %>% dplyr::select(t)  %>% unlist()
  print(Pop2)
  Pop3 <- dataset %>% dplyr::filter(pop==3) %>% dplyr::select(t)  %>% unlist()
  print(Pop3)
  Pop=list(Pop1,Pop2,Pop3)
  return(Pop)
}
Boot <- function(type.dataset,Z.trunc.code,NbSim){
  print(NbSim)
  config=c()
  Data.res <- c()
  Nbconfig <- nrow(Z.trunc.code)
  
  for (c in 1:Nbconfig){
    print(c)
    config.code <- Z.trunc.code[c,]
    # N <- eval(parse(text=paste0("NbSim",type.dataset)))
    if (type.dataset=="Learn"){
      Nb.config.c <- R*0.8
    } else {Nb.config.c <- R*0.2}
    config <- c(config,rep(c,Nb.config.c))
    
    data.c <- purrr::map(List.names.final,~{
      Pop <- c()
      code.names.series <- config.code[which(colnames(Z.trunc.code) %in% .x)]
      eval(parse(text=paste0("Pop=Pop.",.x,".",type.dataset,"[[",code.names.series,"]]")))
      res.t=sample(Pop, Nb.config.c, replace = TRUE)
      return(res.t)
    }) %>% bind_cols() %>% as.data.frame()
    colnames(data.c)=List.names.final
    Data.res = rbind(Data.res,data.c)
  }
  
  Data.res <-Data.res %>% mutate(config=config)
  Data.res$config <- as.factor(Data.res$config)
  return(Data.res)
  
}

PredRule_ET_ErrorTest <- function(DataLearn,DataTest,b,Nbconfig){
  
  config.tot <- 1:Nbconfig
  ##############
  # LDA
  ##############
  modlda<-lda(config~ ., data=DataLearn, CV=FALSE)
  saveRDS(modlda, file = paste0(file_path_Results,'modlda_b',b,significance.level, offset, GE, number.pop,'.rds'))
  
  pred.lda <- predict(modlda,newdata=DataTest,na.action = na.pass)$class 
  err.lda <- mean(pred.lda!=DataTest$config)
  which.misclassif.cluster.lda <- which(pred.lda!=DataTest$config)
  misclassif.cluster.lda <- rbind(DataTest$config[which.misclassif.cluster.lda],pred.lda[which.misclassif.cluster.lda])
  never.pred.lda <- config.tot[!(config.tot %in% pred.lda)]
  
  #####
  #
  cvControl=trainControl(method="cv",number=10)
  #
  ##############
  # Cart 
  ##############
  # Tree.max <- rpart(config ~ ., data = DataLearn,method="class",control=rpart.control(minsplit=2,cp=1e-06))
  # alpha_opt <- Tree.max$cptable[which.min(Tree.max$cptable[,"xerror"]),"CP"]
  # modcart <- prune(Tree.max,cp=alpha_opt)
  # #saveRDS(modcart, file = paste0(file_path_Results,'modcart_b',b,'.rds'))
  # 
  # pred.cart <- predict(modcart,newdata=DataTest,na.action = na.pass,type="class") 
  # err.cart <- mean(pred.cart!=DataTest$config)
  # which.misclassif.cluster.cart <- which(pred.cart!=DataTest$config)
  # misclassif.cluster.cart <- rbind(DataTest$config[which.misclassif.cluster.cart],pred.cart[which.misclassif.cluster.cart])
  # never.pred.cart <- config.tot[!(config.tot %in% pred.cart)]
  # err.cart
  
  modcart = caret::train(config~ ., data = DataLearn, method = "rpart", tuneLength = 10,trControl = cvControl,metric = "Accuracy")
  saveRDS(modcart, file = paste0(file_path_Results,'modcart_b',b,significance.level, offset, GE, number.pop,'.rds'))
  
  pred.cart=predict(modcart, newdata = DataTest)
  err.cart <- mean(pred.cart!=DataTest$config)
  which.misclassif.cluster.cart <- which(pred.cart!=DataTest$config)
  misclassif.cluster.cart <- rbind(DataTest$config[which.misclassif.cluster.cart],pred.cart[which.misclassif.cluster.cart])
  never.pred.cart <- config.tot[!(config.tot %in% pred.cart)]
  err.cart
  
  
  ##############
  # knn
  ##############
  # grid <- c(1,2,3,5,10)
  # error.knn.i <- c()
  # error.knn.i <- purrr::map(grid,~{
  #   pred.knn <- knn(train=DataLearn[,1:5],test=DataTest[,1:5],cl=DataLearn$config,k=.x)
  #   error.knn <- mean(pred.knn!=DataTest$config)
  #   return(error.knn)
  # }) %>% unlist()
  # 
  # k.opt=grid[which.min(error.knn.i)]
  # pred.knn <- knn(train=DataLearn[,1:5],test=DataTest[,1:5],cl=DataLearn$config,k=k.opt)
  modknn = caret::train(config~ ., data = DataLearn, method = "knn", tuneLength = 10,trControl = cvControl,metric = "Accuracy")
  saveRDS(modknn, file = paste0(file_path_Results,'modknn_b',b,significance.level, offset, GE, number.pop,'.rds'))
  
  pred.knn=predict(modknn, newdata = DataTest)
  err.knn <- mean(pred.knn!=DataTest$config)
  which.misclassif.cluster.knn <- which(pred.knn!=DataTest$config)
  misclassif.cluster.knn <- rbind(DataTest$config[which.misclassif.cluster.knn],pred.knn[which.misclassif.cluster.knn])
  never.pred.knn <- config.tot[!(config.tot %in% pred.knn)]
  err.knn
  
  
  
  
  ##############
  # RF
  ##############
  
  # Testé avec le paramètre par défaut pour m=sqrt(p)
  # m=seq(1,4,1)
  # error.rf.i <- c()
  # error.rf.i <- purrr::map(m,~{
  #   foret <- randomForest(config~.,data=DataLearn,mtry=.x,ntree=500)
  #   pred.rf <- predict(foret,newdata=DataTest,na.action = na.pass)
  #   error.rf <- mean(pred.rf!=DataTest$config)
  #   return(error.rf)
  # }) %>% unlist()
  # m.opt=m[which.min(error.rf.i)]
  # modrf <- randomForest(config~.,data=DataLearn,mtry=m.opt)
  # #saveRDS(modrf, file = paste0(file_path_Results,'modrf_b',b,'.rds'))
  # pred.rf <- predict(modrf,newdata=DataTest)
  # err.rf <- mean(pred.rf!=DataTest$config)
  # err.rf
  
  modrf = caret::train(config~ ., data = DataLearn, method = "rf", tuneLength = 4,trControl = cvControl,metric = "Accuracy", importance=TRUE)
  saveRDS(modrf, file = paste0(file_path_Results,'modrf_b',b,significance.level, offset, GE, number.pop,'.rds'))
  pred.rf=predict(modrf, newdata = DataTest)
  err.rf <- mean(pred.rf!=DataTest$config)
  
  which.misclassif.cluster.rf <- which(pred.rf!=DataTest$config)
  misclassif.cluster.rf <- rbind(DataTest$config[which.misclassif.cluster.rf],pred.rf[which.misclassif.cluster.rf])
  never.pred.rf <- config.tot[!(config.tot %in% pred.rf)]
  
  # return
  err.tot <- list(lda=err.lda,cart=err.cart,knn=err.knn,rf=err.rf)
  
  which.misclassif.tot <- list(lda=which.misclassif.cluster.lda,cart=which.misclassif.cluster.cart,knn=which.misclassif.cluster.knn,rf=which.misclassif.cluster.rf)
  
  misclassif.tot <- list(lda=misclassif.cluster.lda,cart=misclassif.cluster.cart,knn=misclassif.cluster.knn,rf=misclassif.cluster.rf)
  
  never.pred.tot <- list(lda=never.pred.lda,cart=never.pred.cart,knn=never.pred.knn,rf=never.pred.rf)
  
  return(list(err.tot=err.tot,
              which.misclassif.tot=which.misclassif.tot,
              misclassif.tot=misclassif.tot,
              never.pred.tot=never.pred.tot,
              modlda=modlda,modcart=modcart,modknn=modknn,modrf=modrf))
  #err.tot=c(err.lda,err.cart,err.knn,err.rf)
  #return(err.tot)
  
}
Boot1 <- function(type.dataset,Z.trunc.code,NbSim){
  config=c()
  Data.res <- c()
  Nbconfig <- nrow(Z.trunc.code)
  
  for (ind.c in 1:Nbconfig){
    print(ind.c)
    config.code <- Z.trunc.code[ind.c,]
    # N <- eval(parse(text=paste0("NbSim",type.dataset)))
    N = NbSim * prob.final[ind.c]
    if (type.dataset=="Learn"){
      Nb.config.c <- N*0.8
    } else {Nb.config.c <- N*0.2}
    config <- c(config,rep(ind.c,Nb.config.c))
    
    data.c <- purrr::map(List.names.final,~{
      Pop <- c()
      code.names.series <- config.code[which(colnames(Z.trunc.code) %in% .x)]
      eval(parse(text=paste0("Pop=Pop.",.x,".",type.dataset,"[[",code.names.series,"]]")))
      res.t=sample(Pop, Nb.config.c, replace = TRUE)
      return(res.t)
    }) %>% bind_cols() %>% as.data.frame()
    colnames(data.c)=List.names.final
    Data.res = rbind(Data.res,data.c)
  }
  
  Data.res <-Data.res %>% mutate(config=config)
  Data.res$config <- as.factor(Data.res$config)
  return(Data.res)
  
}
PredRule_RF<- function(DataLearn,DataTest,b,Nbconfig){
  
  config.tot <- 1:Nbconfig
  
  cvControl=trainControl(method="cv",number=10)
  ##############
  # RF
  ##############
  
  modrf = caret::train(config~ ., data = DataLearn, method = "rf", tuneLength = 4,trControl = cvControl,metric = "Accuracy", importance=TRUE)
  saveRDS(modrf, file = paste0(file_path_Results,'modrf_b',b,significance.level, offset, GE, number.pop,'.rds'))
  pred.rf=predict(modrf, newdata = DataTest)
  err.rf <- mean(pred.rf!=DataTest$config)
  
  which.misclassif.cluster.rf <- which(pred.rf!=DataTest$config)
  misclassif.cluster.rf <- rbind(DataTest$config[which.misclassif.cluster.rf],pred.rf[which.misclassif.cluster.rf])
  never.pred.rf <- config.tot[!(config.tot %in% pred.rf)]
  
  # return
  return(list(err.tot = err.rf,
              which.misclassif.tot = which.misclassif.cluster.rf,
              misclassif.tot = misclassif.cluster.rf,
              never.pred.tot = never.pred.rf,
              modrf = modrf))
  #err.tot=c(err.lda,err.cart,err.knn,err.rf)
  #return(err.tot)
  
}
PredRule_CART<- function(DataLearn,DataTest,b,Nbconfig){
  
  config.tot <- 1:Nbconfig
  
  cvControl=trainControl(method="cv",number=10)
  
  modcart = caret::train(config~ ., data = DataLearn, method = "rpart", tuneLength = 10,trControl = cvControl,metric = "Accuracy")
  saveRDS(modcart, file = paste0(file_path_Results,'modcart_b',b,significance.level, offset, GE, number.pop,'.rds'))
  
  pred.cart=predict(modcart, newdata = DataTest)
  err.cart <- mean(pred.cart!=DataTest$config)
  which.misclassif.cluster.cart <- which(pred.cart!=DataTest$config)
  misclassif.cluster.cart <- rbind(DataTest$config[which.misclassif.cluster.cart],pred.cart[which.misclassif.cluster.cart])
  never.pred.cart <- config.tot[!(config.tot %in% pred.cart)]
  
  # return
  return(list(err.tot = err.cart,
              which.misclassif.tot = which.misclassif.cluster.cart ,
              misclassif.tot = misclassif.cluster.cart,
              never.pred.tot = never.pred.cart,
              mod = modcart))
  #err.tot=c(err.lda,err.cart,err.knn,err.rf)
  #return(err.tot)
  
}
PredRule_knn<- function(DataLearn,DataTest,b,Nbconfig){
  
  config.tot <- 1:Nbconfig
  
  cvControl=trainControl(method="cv",number=10)
  ##############
  # RF
  ##############
  modknn = caret::train(config~ ., data = DataLearn, method = "knn", tuneLength = 10,trControl = cvControl,metric = "Accuracy")
  saveRDS(modknn, file = paste0(file_path_Results,'modknn_b',b,significance.level, offset, GE, number.pop,'.rds'))
  
  pred.knn=predict(modknn, newdata = DataTest)
  err.knn <- mean(pred.knn!=DataTest$config)
  which.misclassif.cluster.knn <- which(pred.knn!=DataTest$config)
  misclassif.cluster.knn <- rbind(DataTest$config[which.misclassif.cluster.knn],pred.knn[which.misclassif.cluster.knn])
  never.pred.knn <- config.tot[!(config.tot %in% pred.knn)]
  err.knn
  
  # return
  return(list(err.tot = err.knn,
              which.misclassif.tot = which.misclassif.cluster.knn,
              misclassif.tot = misclassif.cluster.knn,
              never.pred.tot = never.pred.knn,
              mod = modknn))
  #err.tot=c(err.lda,err.cart,err.knn,err.rf)
  #return(err.tot)
  
}
PredRule_LDA<- function(DataLearn,DataTest,b,Nbconfig){
  
  config.tot <- 1:Nbconfig
  
  # LDA
  ##############
  modlda<-lda(config~ ., data=DataLearn, CV=FALSE)
  saveRDS(modlda, file = paste0(file_path_Results,'modlda_b',b,significance.level, offset, GE, number.pop,'.rds'))
  
  pred.lda <- predict(modlda,newdata=DataTest,na.action = na.pass)$class 
  err.lda <- mean(pred.lda!=DataTest$config)
  which.misclassif.cluster.lda <- which(pred.lda!=DataTest$config)
  misclassif.cluster.lda <- rbind(DataTest$config[which.misclassif.cluster.lda],pred.lda[which.misclassif.cluster.lda])
  never.pred.lda <- config.tot[!(config.tot %in% pred.lda)]
  
  # return
  return(list(err.tot = err.lda,
              which.misclassif.tot = which.misclassif.cluster.lda,
              misclassif.tot = misclassif.cluster.lda,
              never.pred.tot = never.pred.lda,
              mod = modlda))
  #err.tot=c(err.lda,err.cart,err.knn,err.rf)
  #return(err.tot)
  
}