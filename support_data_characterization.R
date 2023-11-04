one.year=365

remove_na_2sides <- function(df, name.series){
  a = which(is.na(df[name.series])== FALSE)
  df = df[c(min(a):(max(a))), ]
  return(df)
}

construct_design <- function(data.df, name.series, break.ind, one.year = 365){
  Data.mod <- data.df %>% dplyr::select(name.series,Date) %>%
    rename(signal=name.series) %>% 
    mutate(complete.time=1:nrow(data.df)) %>% 
    dplyr::select(-Date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one.year),sin",i,"=sin(i*complete.time*(2*pi)/one.year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  n0 = nrow(data.df)
  # Data.mod$right = c(rep(0, break.ind), rep(1, (n0-break.ind)))
  Data.mod$left = rep(1, n0)
  
  return(Data.mod)
}

IGLS <- function(design.m, tol, day.list){
  resi0 = rep(NA, nrow(design.m))
  # call expression
  list.para <- colnames(design.m)[2:dim(design.m)[2]]
  mod.X <-  list.para %>% stringr::str_c(collapse = "+")
  mod.expression <- c("signal","~",mod.X, -1) %>% stringr::str_c(collapse = "")
  # ols
  ols.fit = lm( mod.expression, data = design.m)
  resi0[which(is.na(design.m$signal)==FALSE)] <- ols.fit$residuals
  old.coef = ols.fit$coefficients
  # estimate initial moving variance 
  Y0 = data.frame(date = day.list, residus = resi0)
  w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
  change1 = 10
  i=0
  while (change1 > tol) {
    design.m$w = w0^2
    gls.fit = eval(parse(text=paste0("nlme::gls(",mod.expression,",data=design.m, correlation = NULL, na.action = na.omit, weights=varFixed(value = ~w)",")")))
    change1 = sum((gls.fit$coefficients - old.coef)^2)
    # deg =  cbind(rep(1, nrow(design.m)), as.matrix(design.m[,c(2:9)])) 
    fit.val = as.matrix(design.m[,list.para]) %*% as.matrix(gls.fit$coefficients)
    resi0 = design.m$signal - fit.val
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = 60)
    old.coef = gls.fit$coefficients
    i=1+i
  }
  print(i)
  return(list( coefficients = gls.fit$coefficients, var = w0^2, residual = as.numeric(resi0)))
}

# heteroskedascity
range.var <- function(x, day.list, s){
  df = data.frame(date = day.list, x = x)
  df$y = format(df$date, "%Y")
  df[which(is.na(s)==TRUE),] = NA
  if(all(x==1)){
    NA
  }else{
    anu.min = aggregate(x~y, df, function(z) min(z, na.rm=TRUE))[,2]
    anu.max = aggregate(x~y, df, function(z) max(z, na.rm=TRUE))[,2]
    ifelse(length(anu.max)!=1,  mean( (anu.max-anu.min), na.rm = TRUE), NA)
  }
}

diff.range.var <- function(x, day.list,s){
  df = data.frame(date = day.list, x = x)
  df$y = format(df$date, "%Y")
  df[which(is.na(s)==TRUE),] = NA
  if(all(x==1)){
    NA
  }else{
    anu.min = aggregate(x~y, df, function(z) min(z, na.rm=TRUE))[,2]
    anu.max = aggregate(x~y, df, function(z) max(z, na.rm=TRUE))[,2]
    r = (anu.max- anu.min)
    ifelse(length(anu.max)!=1, (max(r, na.rm = TRUE) - min(r, na.rm = TRUE)), NA)
    
  }
}

p.and.coef <- function(fitARIMA, pq1, nb.or){
  test.sig = lmtest::coeftest(fitARIMA)
  ord = pq1[c(1,3)]
  orde = c(rbind(ord,ord-1))
  orde[which(orde <0)] <- 0
  ind.param = which(orde >0)
  orde[ind.param] <- fitARIMA$coef[1:nb.or]
  p.value <- rep(-1, 4)
  p.value[ ind.param] <- test.sig[,4][1:nb.or]
  return(list(p.value = p.value , coef = orde))
}
# arima 
last_signif <- function(signal, pq, alpha, fit.b){  
  nb.or <- sum(pq)
  pq1 = rep(NA,3)
  while ( identical(as.numeric(pq1), pq) == FALSE) { # iteratively identify the model, stop when the model are the same after the significant check
    pq1 = pq
    if(nb.or==0){
      pandcoef <- list(p.value = rep(-1,4),coef = rep(0,4))
    }else{
      fitARIMA = try(arima( signal, pq, method="ML"), TRUE)
      if (class(fitARIMA) == "try-error"){
        fitARIMA = fit.b
      }
      pandcoef <- p.and.coef(fitARIMA, pq, nb.or)
    }
    pq = check_sig(p.val = pandcoef$p.value, alpha = alpha)
    nb.or <- sum(pq)
  }
  return(list( pq = pq, pandcoef = pandcoef))
}
diff.var <- function(name.test){
  if(name.test == "gps.gps"){
    varname = c("GPS.x", "GPS.y")
  }
  if(name.test == "gps.era"){
    varname = c("GPS.x", "ERAI.x")
  }
  if(name.test == "gps1.era"){
    varname = c("GPS.y", "ERAI.x")
  }
  if(name.test == "gps.era1"){
    varname = c("GPS.x", "ERAI.y")
  }
  if(name.test == "gps1.era1"){
    varname = c("GPS.y", "ERAI.y")
  }
  if(name.test == "era.era"){
    varname = c("ERAI.x", "ERAI.y")
  }
  return(varname)
}

p.and.coef <- function(fitARIMA, pq1, nb.or){
  test.sig = lmtest::coeftest(fitARIMA)
  ord = pq1[c(1,3)]
  orde = c(rbind(ord,ord-1))
  orde[which(orde <0)] <- 0
  ind.param = which(orde >0)
  orde[ind.param] <- fitARIMA$coef[1:nb.or]
  p.value <- rep(-1, 4)
  p.value[ ind.param] <- test.sig[,4][1:nb.or]
  return(list(p.value = p.value , coef = orde))
}
# return significant order
check_sig <- function(p.val, alpha){
  ar.or = length(which(p.val[1:2] >= 0 & p.val[1:2] <= alpha))
  ma.or = length(which(p.val[3:4] >= 0 & p.val[3:4] <= alpha))
  return(c(ar.or, 0, ma.or))
}

model.iden <- function(order){
  model = c()
  if (identical(order, c(1,0,1))){ model = "ARMA(1,1)"}
  else if (identical(order, c(1,0,0))){ model = "AR(1)"}
  else if (identical(order, c(0,0,1))){ model = "MA(1)"}
  else if (identical(order, c(0,0,0))){ model = "White"}
  else if (identical(order, c(2,0,0))){ model = "AR(2)"}
  else if (identical(order, c(2,0,1))){ model = "ARMA(2,1)"}
  else if (identical(order, c(1,0,2))){ model = "ARMA(1,2)"}
  else if (identical(order, c(0,0,2))){ model = "MA(2)"}
  else if (identical(order, c(2,0,2))){ model = "ARMA(2,2)"}
  
  return(model)
}

fit.arma11 <- function(signal.test, significant.level = 0.05){
  fit.b = forecast::auto.arima(signal.test , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean = FALSE,lambda = NULL,
                               max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  
  pq <- arimaorder(fit.b)
  # order.init[k, c((testi*3-2): (testi*3))] <- pq
  options(warn = 2)
  
  refit0 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  pq = refit0$pq
  
  if( any(pq > 1)){
    fit.b = forecast::auto.arima( signal.test, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                  max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
    pq = arimaorder(fit.b)
  }
  
  refit1 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  
  pq = refit1$pq
  return(list(pq = pq, coef = refit1$pandcoef$coef, p = refit1$pandcoef$p.value))
}

extract_param <- function(x){
  if(x =="ARMA(1,1)"){ y = list(rho = 1, theta = 3)}
  if(x == "White"){ y = list(rho = 0, theta = 0)}
  if(x == "AR(1)"){ y = list(rho = 1, theta = 0)}
  if(x == "MA(1)"){ y = list(rho = 0, theta = 3)}
  return(y)
}

convert.name.test <- function(x){
  if (x == "gps.era"){ y = "GPS-ERA"}
  if (x == "gps.era1"){ y = "GPS-ERA'"}
  if (x == "gps1.era"){ y = "GPS'-ERA"}
  if (x == "gps.gps"){ y = "GPS-GPS'"}
  if (x == "era.era"){ y = "ERA-ERA'"}
  if (x == "gps1.era1"){ y = "GPS'-ERA'"}
  return(y)
}

# SLIDING VARIANCE 
RobEstiSlidingVariance.S <- function(Y, name.var, alpha, estimator, length.wind){# require date in the dataset, return std at time t
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = nrow(Y1)
  sigma.est1 = rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    x.idiff = (x.i)
    thre = 30
    if(i < 30|i>(n-30)){thre = 16}
    if(length(na.omit(x.idiff)) <= thre){
      sd <- NA
    }else{
      sd <- my.estimator(estimator = estimator, x.idiff)
    }
    sigma.est1[i] <- sd
  }
  # linear regression of the variance for gaps  MAYBE REPLACE BY INTERPOLATION FUNCTION
  s = sigma.est1
  if (sum(is.na(s)) != 0 & sum(is.na(s)) != length(s)){
    ts.s = c(1:n)
    na.ind = which(is.na(s))
    if(na.ind[1] == 1){
      ind.stop = which(is.na(s)==FALSE)[1]-1
      na.ind <- na.ind[-c(1:ind.stop)]
    }else if (is.na(s[n]) == 1){
      m = which(is.na(s)==FALSE)
      ind.start = m[length(m)]
      na.ind <- na.ind[-which(na.ind %in% c(ind.start:n))]
    }
    s[na.ind] <- approx(ts.s, s, xout=na.ind)$y
  }
  sigma.est = s[which(Y1$date %in% Y$date)]
  return(sigma.est)
}
my.estimator <- function(estimator,x){
  x1 = na.omit(x)
  if(estimator == "mad"){
    n0 = length(x1)
    f <- qnorm(3/4)*(1-0.7612/n0 - 1.123/(n0^2))
    y = mad(x1, constant = 1/f)
  }else if(estimator == "Qn"){
    y = robustbase::Qn(x1)
  }else if(estimator == "Sca"){
    y = robustbase::scaleTau2(x1, consistency = "finiteSample")
  }
  return(y)
}

fit.arima <- function(signal.test, significant.level = 0.05){
  fit.b = forecast::auto.arima(signal.test , d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                               max.p = 2, max.q = 2, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
  
  pq <- forecast::arimaorder(fit.b)
  # order.init[k, c((testi*3-2): (testi*3))] <- pq
  # options(warn = 2)
  
  refit0 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  pq = refit0$pq
  
  if( any(pq > 1)){
    fit.b = forecast::auto.arima( signal.test, d = 0, ic = "bic", seasonal = FALSE, stationary = TRUE, allowmean =FALSE,lambda = NULL,
                                  max.p = 1, max.q = 1, start.p = 0, trace = FALSE, allowdrift = FALSE,  approximation=FALSE)
    pq = forecast::arimaorder(fit.b)
  }
  test.pq = pq
  refit1 = last_signif(signal = signal.test, pq, alpha = significant.level, fit.b = fit.b)
  
  pq = refit1$pq
  if(identical(as.numeric(test.pq),pq) == FALSE){print(c(test.pq,pq))}
  return(list(pq = pq, coef = refit1$pandcoef$coef, p = refit1$pandcoef$p.value))
}

plot_for_papar <- function(arma_order, arma_coef, list_infor, name_six_diff){
  order_arma = read.table(file = paste0(path_results, "order_arma.txt"), header = TRUE)
  coef_arma = read.table(file = paste0(path_results, "coef_arma.txt"), header = TRUE)
  
  list_model = c("White", "AR(1)", "MA(1)", "ARMA(1,1)")
  
  text1 = "Distance < 50 km"
  text2 = "Distance > 50 km"
  
  length_data = nrow(order_arma)
  six_model = data.frame(matrix(NA, ncol = 6, nrow = length_data))
  for (i in 1:6) {
    six_model[,i] = sapply(c(1:length_data), function(x) {
      model.iden(as.numeric(unlist(order_arma[x, (3*i-2):(3*i)])))
    })
  }
  colnames(six_model) <- name_six_diff
  
  six_values = c()
  for (i in 1:6) {
    value_count = sapply(c(list_model), function(x) length(which(six_model[,i] == x)))
    six_values <- c( six_values, value_count)
  }
  
  all_model <- cbind(list_infor[,c(1:3)], six_model) %>%
    mutate(main_nb = if_else(main < nearby, 
                             paste(main, nearby, sep = ""), 
                             paste(nearby, main, sep = ""))) 
  
  consistency_check <- cbind(list_infor[,c(1:3)], six_model) %>%
    mutate(main_nb = if_else(main < nearby, 
                             paste(main, nearby, sep = ""), 
                             paste(nearby, main, sep = "")))  %>%
    group_by(main_nb) %>%
    summarise(Consistent = all(GPS_ERA == GPS_ERA[1]) &
                all(GPS_GPS1 == GPS_GPS1[1]) &
                all(GPS_ERA1 == GPS_ERA1[1]) &
                all(ERA_ERA1 == ERA_ERA1[1]) &
                all(GPS1_ERA1 == GPS1_ERA1[1]) &
                all(GPS1_ERA == GPS1_ERA[1]),
              Count = n()
    ) %>%
    ungroup()
  
  # Function to check consistency for each model
  model_inconsistencies <- all_model %>%
    group_by(main_nb) %>%
    summarise(
      Inconsistency_GPS_ERA = n_distinct(GPS_ERA),
      Inconsistency_GPS_GPS1 = n_distinct(GPS_GPS1),
      Inconsistency_GPS_ERA1 = n_distinct(GPS_ERA1),
      Inconsistency_ERA_ERA1 = n_distinct(ERA_ERA1),
      Inconsistency_GPS1_ERA1 = n_distinct(GPS1_ERA1),
      Inconsistency_GPS1_ERA = n_distinct(GPS1_ERA)
    )
  
  res_plot = data.frame(series = rep(list_name_test, each = 4), 
                        mod = rep(list_model, 6),
                        value = six_values, 
                        n = rep(length_data,24))
  res_plot$pct = res_plot$value/res_plot$n*100
  res_plot$series = factor(res_plot$series, 
                           levels=reoder_list_name)
  res_plot$mod = factor(res_plot$mod, 
                        levels=list_model)
  
  p1 <- ggplot(res_plot, aes(fill=mod, y=pct, x=series, label = value)) + 
    geom_bar(position="dodge", stat="identity", width = 0.5)+theme_bw()+ 
    xlab("") + ylab("Percentage")+
    labs(tag = "(a)", subtitle = text1) +
    geom_text(position = position_dodge(width = .5),    # move to center of bars
              vjust = -0.5,    # nudge above top of bar
              size = 1)+
    ylim(c(0,100))+
    theme(axis.text.x = element_text(size = 4.5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
          axis.title = element_text(size = 5), legend.key.size = unit(0.2, "cm"), 
          plot.tag = element_text(size = 6) , plot.subtitle = element_text(size = 6),
          legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))
  
  ggsave(paste0(path_results,"attribution/Datacharacterization_autoarima_test.jpg" ), plot = p1, width = 14, height = 5, units = "cm", dpi = 1200)
  
  
  
}


