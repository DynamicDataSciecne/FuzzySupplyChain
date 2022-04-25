# Auxiliary function for Fuzzy SC paper
# Author: Erfanul Hoque
# Date: July 2021


######################################################################
# Auxillary functions
######################################################################
#library(tseries)

######################################################################
# Generate AR(1) model # D_t + \mu = \phi (D_{t-1} - \mu) + a_t
######################################################################
# Simulate AR(p) model
AR_sim <- function(t=100, mean = 0, ar.coef = 0.75, ma.coef = NULL, order = NULL, 
                   error = c("norm", "t"), sd = 1, df=3.5){
  error <- match.arg(error)
  # simulating data from AR(p) model
  D.t <- if (error == "norm") {
    # Note: arima.sim() will assume normally distributed errors with SD = 1
    arima.sim(n = t, list(order = order, ar = ar.coef, ma=ma.coef), sd = sd) + mean
  } else { #long-tailed
    arima.sim(n = t, list(order = order, ar = ar.coef, ma = ma.coef), rand.gen = function(t, ...) 
      rt(t, df = df)) + mean
  }
  return(D.t)
}

# data generation from AR(1)
#set.seed(2021)
#D.t <- AR_sim(t=400, mean = 500, ar.coef = 0.8, order = c(1, 0, 0), error = "norm", sd= 3) # AR(1)
#D.t <- AR_sim(t=400, mean = 500, ar.coef = 0.8, order = c(1, 0, 0), error = "t", df=3.5) # AR(1)
# head(D.t)
## Check variance for normal
# sd(D.t)
# sqrt(3^2/(1-0.8^2))

######################################################################
# Estimate AR(p) model : D_t + \mu = \phi (D_{t-1} - \mu) + a_t
######################################################################
# fit the AR model to simulated data using "arima" function 
# arima: ARIMA Modelling of Time Series

AR_est <- function(data, order = c(1, 0, 0), include.mean = TRUE, # for "arima" function
                   #include.constant = FALSE, # for "Arima" function
                   method = c("CSS-ML", "ML", "CSS")){
  # method: fitting method- maximum likelihood or minimize conditional sum-of-squares. 
  # The default (unless there are missing values) is to use conditional-sum-of-squares 
  # to find starting values, then maximum likelihood.
  res <- arima(x = data, order = order, include.mean = TRUE)
  #res <- forecast::Arima(D.t, order = c(1, 0, 0), include.constant = TRUE)
  return(res)
}

######################################################################
# L-step ahead fuzzy forecast for simulation
#####################################################################
Fuzzy_forecastARsim <- function(data, L=5, alpha.cut) {
  
  #Fit data and fuzzifying parameters
  res <- AR_est(data = data, order = c(1, 0, 0))
  #Extract parameters
  Mu <- as.numeric(res$coef)[2]
  phi <- as.numeric(res$coef)[1]
  # L-step ahead MMSE forecast for AR1
  errorvar <- res$sigma2
  #MMSE.hat <- phi*data[l]
  # calculate sum of L-step ahead forecast
  d.t <- data[length(data)]
  MMSE.hat <- (((phi/(1-phi))*(1-phi^L))*d.t)
  var_mmse <- ((phi^2/(1-phi)^2)*(1-phi^L)^2)* var(data) # var(data): sample var of original demand
  
  MMSE_fore.L <- sapply(1:length(alpha.cut), function(i) (MMSE.hat + qnorm(alpha.cut[i]/2, 0, 1) * sqrt(var_mmse)))
  MMSE_fore.U <- sapply(1:length(alpha.cut), function(i) (MMSE.hat - qnorm(alpha.cut[i]/2, 0, 1) * sqrt(var_mmse)))
  
  #MMSE_fore.L <- sapply(1:length(alpha.cut), function(i) (MMSE.hat + qnorm(alpha.cut[i]/2, 0, 1) * sqrt(errorvar)))
  #MMSE_fore.U <- sapply(1:length(alpha.cut), function(i) (MMSE.hat - qnorm(alpha.cut[i]/2, 0, 1) * sqrt(errorvar)))
    
  results <- data.frame(alpha.cut, MMSE_fore.L, MMSE_fore.U)
  return(results)
}

######################################################################
# L-step ahead fuzzy forecast for data
#####################################################################
Fuzzy_forecastARdata <- function(data, L=5, alpha.cut) {
  
  #Fit data and fuzzifying parameters
  res <- AR_est(data = data, order = c(1, 0, 0))
  #Extract parameters
  Mu <- as.numeric(res$coef)[2]
  phi <- as.numeric(res$coef)[1]
  # L-step ahead MMSE forecast for AR1
  errorvar <- res$sigma2
  #MMSE.hat <- phi*data[l]
  # calculate sum of L-step ahead forecast
  d.t <- data[length(data)]
  MMSE.hat <- (((phi/(1-phi))*(1-phi^L))*d.t)
  var_mmse <- ((phi^2/(1-phi)^2)*(1-phi^L)^2)* var(data) # var(data): sample var of original demand
  
  #MMSE_fore.L <- sapply(1:length(alpha.cut), function(i) (MMSE.hat + qnorm(alpha.cut[i]/2, 0, 1) * sqrt(var_mmse)))
  #MMSE_fore.U <- sapply(1:length(alpha.cut), function(i) (MMSE.hat - qnorm(alpha.cut[i]/2, 0, 1) * sqrt(var_mmse)))
  
  MMSE_fore.L <- sapply(1:length(alpha.cut), function(i) (MMSE.hat + qnorm(alpha.cut[i]/2, 0, 1) * sqrt(errorvar)))
  MMSE_fore.U <- sapply(1:length(alpha.cut), function(i) (MMSE.hat - qnorm(alpha.cut[i]/2, 0, 1) * sqrt(errorvar)))
  
  results <- data.frame(alpha.cut, MMSE_fore.L, MMSE_fore.U)
  return(results)
}

########################## DDWMA method ###############################

#######################################################################
# Function to calculate lag L (L= 1, 2, 3, 4, 5,...) values
#######################################################################
lag_data <- function(data, L){
  data <- as.vector(data)
  lag.data <- sapply(1:L, function(l) dplyr::lag(data, l))
  return(lag.data)
}
# Calculate lag values
# lag.data <- lag_data(data = data, L = L)

#######################################################################
# Function to calculate the optimal weights over number of simulations
######################################################################

ML_approach <- function(data, L = 5, replication = 500){
  L <- L
  # Calculate lag values
  lag.data <- lag_data(data = data, L = L)
  # Do the simulation for number of replication times
  # Creating a matrix to store the weights
  all_wts <- matrix(nrow = replication, ncol = ncol(lag.data))
  # Creating an empty vector to store demand forecast
  dem_fore <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
  # Creating an empty vector to store Standard deviation
  meanFore <- vector('numeric', length = replication)
  riskFore <- vector('numeric', length = replication)
  # Creating an empty vector to store Sharpe Ratio
  sharpe_ratio <- vector('numeric', length = replication)
  
  error <- matrix(0, nrow = nrow(lag.data) - L, ncol = replication)
  # Creating an empty vector to store results
  SSE <- vector('numeric', length = replication)
  #all_wts <- vector('numeric', length = replication)
  
  # Do the simulation
  for (i in 1:replication) {
    set.seed(2021 + i)
    # Create random weights first.
    wts <- runif(n = L)
    wts <- wts/sum(wts) # standardized the weight to get sum equal 1
    all_wts[i,] <- wts
    # demand forecast (forecast is valid only when t >= L+1)
    demand_fore <- sapply((L+1):nrow(lag.data), function(l) sum(wts * lag.data[l,]))
    # Storing demand forecast values
    dem_fore[,i] <- demand_fore
    # Calculate the risk
    meanForeF <- mean(demand_fore)
    meanFore[i] <- meanForeF 
    dem_risk <- sd(demand_fore)
    riskFore[i] <- dem_risk
    sharpe_ratio[i] <- meanForeF / dem_risk
    
    error[,i] <- (data[-c(1:L)] - dem_fore[,i])
    SSE[i] <- sum(error[,i]^2) 
    
  }
  # Storing the values in the table
  res1 <- data.frame(Weight = all_wts, SSE = SSE)
  res_values <- list(results = res1, demandForecast = dem_fore, error= error, 
                     MeanDemand = meanFore, Risk = riskFore, SRatio = sharpe_ratio,
                     Lagvalues = lag.data) 
  return(res_values)
}

######################################################################
# L-step ahead fuzzy forecast of DDWMA
#####################################################################
# Calculate L-step ahead fuzzy forecast
# D_{t+1}_bar(L) = DDWMA + z_\alpha/2 * se(DDWMA)

Fuzzy_fore.DDWMA <- function(data, L=5, alpha.cut) {
  
  #Fit data and fuzzifying parameters
  res <- AR_est(data = data, order = c(1, 0, 0))
  #Extract parameters
  Mu <- as.numeric(res$coef)[2]
  phi <- as.numeric(res$coef)[1]
  
  all_res <- ML_approach(data = data, L = L, replication = 5000)
  Optm.res <- all_res$results
  min_sd <-  Optm.res[which.min(Optm.res$SSE),]
  L <- L
  Optweights <- min_sd[1:L]
  Optlagvalues <- all_res$Lagvalues[length(data),]
  phi.fit <- phi
  minMSE <- min_sd$SSE/length(data)
  minvar <- all_res$Risk[which.min(Optm.res$SSE)]^2
  
  # calculate sum of L-step ahead forecast  
  Sfore_ML <- sum(Optweights * Optlagvalues)
  
  # L-step ahead fuzzy DDWMA forecast
  #DDWMA_fore.L <- sapply(1:length(alpha.cut), function(i) (Sfore_ML + qnorm(alpha.cut[i]/2, 0, 1) * sqrt(minMSE)))
  #DDWMA_fore.U <- sapply(1:length(alpha.cut), function(i) (Sfore_ML - qnorm(alpha.cut[i]/2, 0, 1) * sqrt(minMSE)))

  DDWMA_fore.L <- sapply(1:length(alpha.cut), function(i) (Sfore_ML + qnorm(alpha.cut[i]/2, 0, 1) * sqrt(minvar)))
  DDWMA_fore.U <- sapply(1:length(alpha.cut), function(i) (Sfore_ML - qnorm(alpha.cut[i]/2, 0, 1) * sqrt(minvar)))

  results <- data.frame(alpha.cut, DDWMA_fore.L, DDWMA_fore.U)
  return(results)
}


######################################################################
# L-step ahead fuzzy Bollinger Bands forecast using SMA
#####################################################################

fuzzyBB_sma <- function(data, L = 5, dist = c("norm", "t"), alpha.cut){
  tt <- length(data)
  sma.values <- SMA(data, n=L)
  ss <- sma.values
  sd.sma <- rollapply(data, width = L, FUN = sd, by.column = TRUE, fill = NA, align = "right")
  
  # Bollinger band using 2-sigma
  bb.lower.sd <- sma.values - 2*sd.sma
  bb.upper.sd <- sma.values + 2*sd.sma
  
  #calculate residual of demand
  #ss <- c(rep(NA,5), sma.values) 
  res <- na.omit(data - ss)
  # Sign correlation
  rho.cal <- function(y) cor(y-mean(y), sign(y-mean(y)))
  rho <- rho.cal(res)
  
  if(dist == "t"){
    nu.fun <- function (x) rho*(x-1)*beta(x/2,1/2)-2*sqrt(x-2)
    nu <- uniroot(nu.fun, c(2, 15))$root
  }else {nu <- 4}
  
  #data driven volatility estimate (DDVE)
  vol.cal <-function(y, rho){return(mean(abs(y-mean(y)))/rho)}
  vol.sma <- rollapply(data, width = L, FUN = vol.cal, rho = rho, by.column = TRUE, fill = NA,align = "right")
 # vol.sma <- sapply((L+1):length(data), function(t) mean(abs(data[(t-5+1):t]- ss[t]))/rho)
  
  # L-step ahead fuzzy SMA forecast
  if(dist == "norm"){
    bbt.lower.vol1 <- sapply(1:length(alpha.cut), function(i) 
      (sma.values[tt] + qnorm(alpha.cut[i]/2, 0, 1)*vol.sma[tt]))
    bbt.upper.vol1 <- sapply(1:length(alpha.cut), function(i) 
      (sma.values[tt] - qnorm(alpha.cut[i]/2, 0, 1)*vol.sma[tt]))
  } else{
    bbt.lower.vol1 <- sapply(1:length(alpha.cut), function(i)
      (sma.values[tt] - qstd((1 - (alpha.cut[i]/2)), mean = 0, sd = 1, nu=nu)*vol.sma[tt]))
    bbt.upper.vol1 <- sapply(1:length(alpha.cut), function(i)
      (sma.values[tt] + qstd((1 - (alpha.cut[i]/2)), mean = 0, sd = 1, nu=nu)*vol.sma[tt]))
  }
  
  BB.res <- as.data.frame(cbind(bbt.lower.vol1, bbt.upper.vol1))
  colnames(BB.res)<-c("BBL_SMA", "BBU_SMA")
  
  return(BB.res)
}  


#######################################################################
# Function to calculate the optimal weights for volatility forecast 
######################################################################
# data: the original data
ML_approach_Vol <- function(data, L = 5, replication = 5000){
  L <- L
  # First calculate the optimal weighted demand forecast based on DDWMA
  all_res <- ML_approach(data = data, L = 5, replication = 5000)
  DDWMA.values <- all_res$demandForecast[,which.min(all_res$results$SSE)]
  opt.weight.demand <- all_res$results[which.min(all_res$results$SSE),]
  opt.error.demand <- all_res$error[,which.min(all_res$results$SSE)]
  # Extract the residuals of demand for volatility estimate
  res <- all_res$error[,which.min(all_res$results$SSE)]
  # Sign correlation and d.f
  rho.cal <- function(y) cor(y-mean(y), sign(y-mean(y)))
  rho <- rho.cal(res)
  
  nu.fun <- function (x) rho*(x-1)*beta(x/2,1/2)-2*sqrt(x-2)
  nu <- uniroot(nu.fun, c(2, 15))$root
  
  # convert residual to standardized using rho to get original volatility
  #res.original <- res/rho 
  res.original.vol <- abs(res)/rho 
  
  # Calculate L-lag values for volatility forecast (L window size)
  lag.data.res <- lag_data(data = res.original.vol, L = L)
  
  # Creating a matrix to store the weights fro volatility
  all_wts <- matrix(nrow = replication, ncol = L)
  # Creating an empty vector to store volatility forecast
  vol_fore <- matrix(0, nrow = nrow(lag.data.res) - L, ncol = replication)
  # Creating an empty vector to store Standard deviation
  mean.vol <- vector('numeric', length = replication)
  sd.vol <- vector('numeric', length = replication)
  # Creating an empty vector to store FESS
  error.vol <- matrix(0, nrow = nrow(lag.data.res) - L, ncol = replication)
  #Creating an empty vector to store SSE
  SSE.vol <- vector('numeric', length = replication)
  
  # Do the simulation to get optimal weight for volatility based on standardized residual data
  for (i in 1:replication) {
    set.seed(2021 + i)
    # Create random weights first.
    wts <- runif(n = L)
    wts <- wts/sum(wts) # standardized the weight to get sum equal 1
    all_wts[i,] <- wts
    # calculate volatility forecast based on L window size (note: first L observations will be missing)
    vol.fore <- sapply((L+1):nrow(lag.data.res), function(t) sum(wts * lag.data.res[t,]))
    # Storing demand forecast values
    vol_fore[,i] <- vol.fore
    
    # Calculate the risk
    meanVolF <- mean(vol.fore)
    mean.vol[i] <- meanVolF 
    vol.risk <- sd(vol.fore)
    sd.vol[i] <- vol.risk
    
    # Calculate one-step ahead forecast error (FESS) to choose optimal weight
    error.vol[,i] <- (res.original.vol[-c(1:L)] - vol_fore[,i])
    SSE.vol[i] <- sum(error.vol[,i]^2) 
  }
  
  # Storing the values in the data table
  res1 <- data.frame(Weight = all_wts, SSE = SSE.vol)
  res_values <- list(results = res1, VolForecast = vol_fore, error.vol = error.vol, 
                     demandForecast=DDWMA.values, opt.weight.demand = opt.weight.demand,
                     error.demand = opt.error.demand,
                     meanVolatility = mean.vol, sdVolatility = sd.vol, est.df=nu, rho=rho) 
  return(res_values)
}



######################################################################
# L-step ahead fuzzy Bollinger Bands forecast using DDWMA
#####################################################################

fuzzyBB_DDWMA <- function(data, L = 5, dist = c("norm", "t"), alpha.cut){
  tt <- length(data)
  ###
  all_res <- ML_approach_Vol(data = data, L = 5, replication = 5000)
  # Obtain the optimal weight based on minimum FESS
  Optm.res <- all_res$results
  min_vol <- Optm.res[which.min(Optm.res$SSE),]  # subtract the results to calculate the optimal weight
  # optimal volatility and weights
  vol.for <- all_res$VolForecast[,which.min(Optm.res$SSE)]
  opt.weight.vol <- min_vol
  error.vol <- all_res$error.vol[,which.min(Optm.res$SSE)]
  # optimal demand forecast and optimal weight 
  DDWMA.values <- all_res$demandForecast[-c(1:L)]
  opt.weight.demand <- all_res$opt.weight.demand
  error.demand <- all_res$error.demand[-c(1:L)]
  # Sign correlation
  rho <- all_res$rho
  nu <- all_res$est.df
  
  # L-step ahead fuzzy DDWMA forecast
  if(dist == "norm"){
    bbt.lower.vol1 <- sapply(1:length(alpha.cut), function(i) 
      (DDWMA.values[tt-2*L] + qnorm(alpha.cut[i]/2, 0, 1)*vol.for[tt-2*L]))
    bbt.upper.vol1 <- sapply(1:length(alpha.cut), function(i) 
      (DDWMA.values[tt-2*L] - qnorm(alpha.cut[i]/2, 0, 1)*vol.for[tt-2*L]))
  } else{
    bbt.lower.vol1 <- sapply(1:length(alpha.cut), function(i)
      (DDWMA.values[tt-2*L] - qstd((1 - (alpha.cut[i]/2)), mean = 0, sd = 1, nu=nu)*vol.for[tt-2*L]))
    bbt.upper.vol1 <- sapply(1:length(alpha.cut), function(i)
      (DDWMA.values[tt-2*L] + qstd((1 - (alpha.cut[i]/2)), mean = 0, sd = 1, nu=nu)*vol.for[tt-2*L]))
  }
  
  BB.res <- as.data.frame(cbind(bbt.lower.vol1, bbt.upper.vol1))
  colnames(BB.res)<-c("BBL_DDWMA", "BBU_DDWMA")
  
  return(BB.res)
}  

