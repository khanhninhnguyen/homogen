#' FGLS  NOTE THAT THIS IS ONLY USED FOR FGLS
#' SHOULD NOT BE LOADED WHEN DOING CHARACTERIZATION OR OTHERS 

library(nlme)
# Function are used 
# construct design matrix -- need to be full date
construct_design <- function(data_df, name_series, break_ind, one_year){
  data_df <- tidyr::complete(data_df, Date = seq(min(data_df$Date), max(data_df$Date), by = "day"))
  Data.mod <- data_df %>% dplyr::select(name_series,Date) %>%
    rename(signal=name_series) %>% 
    mutate(complete.time=1:nrow(data_df)) %>% 
    dplyr::select(-Date)
  for (i in 1:4){
    eval(parse(text=paste0("Data.mod <- Data.mod %>% mutate(cos",i,"=cos(i*complete.time*(2*pi)/one_year),sin",i,"=sin(i*complete.time*(2*pi)/one_year))")))
  }
  Data.mod <- Data.mod %>% dplyr::select(-complete.time)
  n0 = nrow(data_df)
  Data.mod$right = c(rep(0, break_ind), rep(1, (n0-break_ind)))
  Data.mod$left = rep(1, n0)
  
  return(Data.mod)
}
# construct the correlation matrix 
cor_matrix <- function(phi, theta, n.full, n.true){
  if(phi ==0 & theta ==0){
    y = rep(1, n.full)
  }else{
    y = ARMAacf(ar = phi, ma = theta, lag.max = n.true-1, pacf = FALSE)
  }
  return(toeplitz(y))
}

GLS <- function(phi, theta, var.t, design.matrix){
  # compute the variance covariance matrix 
  ind1 = which(is.na(design.matrix$signal)==FALSE)
  var.t.na = var.t[ind1]
  var.matrix = diag(sqrt(var.t.na))
  cor.matrix = cor_matrix(phi, theta, n.full = nrow(design.matrix), n.true = nrow(na.omit(design.matrix)))
  if(phi==0&theta==0){
    cov.var= diag(var.t.na)
  }else{
    cov.var0 = var.matrix  %*%  cor.matrix 
    cov.var = cov.var0 %*%  var.matrix
  }
  
  # estimate
  X = as.matrix(design.matrix%>% dplyr::select(-signal))
  X = X[ind1,]
  term2 = solve(cov.var)
  term1 = t(X) %*% term2 %*% X
  term3 = solve(term1)
  beta = term3 %*% t(X) %*% term2  %*% (as.matrix(design.matrix$signal[ind1]))
  var.beta = term3 
  
  # form the frame of result as in the ols 
  residual = design.matrix$signal[ind1] - X%*%beta
  t.val = beta/sqrt((diag(var.beta)))
  p.val = round(pnorm(-abs(t.val), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)
  fit.gls = data.frame(Estimate = beta)
  fitted.gls = X%*%beta
  fit.gls$`Std. Error` = sqrt(diag(var.beta))
  fit.gls$`t value` = t.val
  fit.gls$`Pr(>|t|)` = p.val
  
  return(list(Coefficients = beta, t.table = fit.gls, vcov = var.beta, residual = residual, fitted = fitted.gls))
}

chooseparam <- function(noise.model, arima.fit){
  if(noise.model[1]==1 & noise.model[3]==1){
    phi = arima.fit$coef[1]
    theta = arima.fit$coef[2]
  }else if(noise.model[1]==1 & noise.model[3]==0){
    phi = arima.fit$coef[1]
    theta = 0
  }else if(noise.model[1]==0 & noise.model[3]==1){
    phi = 0
    theta = arima.fit$coef[1]
  }else if(noise.model[1]==0 & noise.model[3]==0){
    phi = 0
    theta = 0
  }
  return(list(phi = phi, theta = theta))
}

FGLS1 <- function(design.m, tol, day.list, noise.model, length.wind0){
  start_time <- Sys.time()
  resi0 = rep(NA, nrow(design.m))
  ind1 = which(is.na(design.m$signal)==FALSE)
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
  w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = length.wind0)
  change1 = 10
  change2 = 10
  change3 =10
  coef.old= list(phi = 0, theta =0)
  i=0
  # print(i)
  while (change3 > tol) {
    # design.m$w0 = w0
    # estimate phi and theta 
    # fit.wls = WLS(var.t = w0, design.matrix = design.m)
    design.m$signal[which(is.na(w0)==TRUE)] = NA
    design.mWLS = design.m
    design.mWLS$w0 = w0
    fit.wls = lm(mod.expression, data = design.mWLS, weights=w0)
    # print(i*10)
    resi0 <- design.m$signal -as.matrix(design.m[,-1]) %*% as.matrix(fit.wls$coefficients)
    normalized.res = resi0/sqrt(w0)
    # noise.model = fit.arima(normalized.res )
    arima.fit = arima(x =normalized.res, order = noise.model, include.mean = FALSE)
    # # GLS
    coef.arma = chooseparam(noise.model = noise.model, arima.fit = arima.fit)
    fit.gls = GLS(phi = coef.arma$phi, theta = coef.arma$theta, var.t = w0, design.matrix = design.m)
    # use th gls function 
    # fit.gls <- eval(parse(text=paste0("gls(",mod.expression,",data=design.m,correlation =", cor.struct, "na.action=na.omit,weights=varFixed(value = ~w0)",")")))
    # update moving variance 
    fit.val = fit.gls$fitted
    resi0 = rep(NA, nrow(design.m))
    ind1 = which(is.na(design.m$signal)==FALSE)
    resi0[ind1] = fit.gls$residual
    t.table = fit.gls$t.table
    Y0 = data.frame(date = day.list, residus = resi0)
    w0 = RobEstiSlidingVariance.S(Y = Y0, name.var = "residus", alpha = 0, estimator = "Sca", length.wind = length.wind0)
    normal.beta = (fit.gls$Coefficients - old.coef)/(sqrt(diag(fit.gls$vcov)))
    change2 = change1-sum((normal.beta)^2)
    change1 = sum((normal.beta)^2)
    # change3
    old.coef = fit.gls$Coefficients
    i=1+i
    j=0
    if(i>10){
      if(change2<(tol/100)){j=1}
      break
    }
    if(identical(noise.model,c(0,0,1))){
      change3 = abs(coef.arma$theta - coef.old$theta)
    }else{
      change3 = abs(coef.arma$phi - coef.old$phi)
    }
    coef.old = coef.arma                
  }
  print(i)
  end_time <- Sys.time()
  design.m$residual = resi0
  design.m$norm.res = normalized.res
  design.m$date = day.list
  return(list( coefficients = fit.gls$Coefficients, var = w0, residual = resi0, fit = fit.val, t.table=t.table, coef.arma = coef.arma,  
               i=i, j = j, change1= change1, all.out = fit.gls, t = (end_time - start_time), design.matrix = design.m))
}
remove_na_2sides <- function(df, name.series){
  a = which(is.na(df[name.series])== FALSE)
  df = df[c(min(a):(max(a))), ]
  return(df)
}
