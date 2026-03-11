############################################################
# Simulation Study for VREHT Detection Function
# True detection model: 0.5 * HN(20) + 0.5 * HR(10,7)
# Competing models: Half-Normal (HN), Hazard-Rate (HR), VREHT
# Author: Gajanan Patil
############################################################

rm(list=ls())

############################################
# SETTINGS
############################################

w <- 50                 # truncation distance

sigma_HN_true <- 20     # true scale parameter for HN
sigma_HR_true <- 10     # true scale parameter for HR
b_true <- 7             # HR shape parameter
pi_true <- 0.5          # mixture proportion

nsim <- 500             # number of simulations
n <- 100                # sample size
eps <- 1e-5             # numerical derivative step

breaks <- seq(0,w,length.out=6)  # grouped intervals

############################################
# DETECTION FUNCTIONS
############################################

# Half-Normal detection function
g_HN <- function(x,sigma){
  exp(-(x^2)/(2*sigma^2))
}

# Hazard-Rate detection function
g_HR <- function(x,sigma,b){
  1-exp(-(x/sigma)^(-b))
}

# VREHT detection function
g_VREHT <- function(x,sigma,a){
  1-(tanh(x/sigma))^a
}

############################################
# TRUE MIXTURE MODEL
############################################

# mixture detection function used to generate data
g_mix <- function(x){
  pi_true*g_HN(x,sigma_HN_true) +
  (1-pi_true)*g_HR(x,sigma_HR_true,b_true)
}

# true average detection probability
true_pa <- integrate(function(x) g_mix(x),1e-6,w)$value/w

############################################
# SIMULATE DISTANCES
############################################

simulate_mix <- function(n){

  y <- numeric()

  while(length(y)<n){

    x <- runif(n,0,w)
    u <- runif(n)

    accept <- u < g_mix(x)

    y <- c(y,x[accept])
  }

  y[1:n]
}

############################################
# INTERVAL PROBABILITY
############################################

interval_prob <- function(gfun,param){

  k <- length(breaks)-1
  p <- numeric(k)

  denom <- integrate(function(x) gfun(x,param),1e-6,w)$value

  for(i in 1:k){

    num <- integrate(function(x) gfun(x,param),
                     max(breaks[i],1e-6),
                     breaks[i+1])$value

    p[i] <- num/denom
  }

  p
}

############################################
# STORAGE
############################################

results <- data.frame()
valid_sim <- 0
sim_attempt <- 0

############################################
# SIMULATION LOOP
############################################

while(valid_sim < nsim){

  sim_attempt <- sim_attempt + 1

  y <- simulate_mix(n)

  grp <- cut(y,breaks,include.lowest=TRUE)
  counts <- as.numeric(table(grp))

############################################
# HALF-NORMAL MODEL
############################################

  loglik_HN <- function(par){

    sigma <- par
    p <- interval_prob(function(x,p) g_HN(x,p),sigma)

    -sum(counts*log(p))
  }

  fit_HN <- try(optim(20,loglik_HN,hessian=TRUE),silent=TRUE)
  if(inherits(fit_HN,"try-error")) next

  sigma_HN <- fit_HN$par

  cov_HN <- try(solve(fit_HN$hessian),silent=TRUE)
  if(inherits(cov_HN,"try-error")) next

  pa_HN <- integrate(function(x) g_HN(x,sigma_HN),1e-6,w)$value/w

############################################
# HAZARD-RATE MODEL
############################################

  loglik_HR <- function(par){

    sigma <- par[1]
    b <- par[2]

    p <- interval_prob(function(x,p) g_HR(x,p[1],p[2]),c(sigma,b))

    -sum(counts*log(p))
  }

  fit_HR <- try(optim(c(20,5),loglik_HR,hessian=TRUE),silent=TRUE)
  if(inherits(fit_HR,"try-error")) next

  sigma_HR <- fit_HR$par[1]
  b_HR <- fit_HR$par[2]

  cov_HR <- try(solve(fit_HR$hessian),silent=TRUE)
  if(inherits(cov_HR,"try-error")) next

  pa_HR <- integrate(function(x) g_HR(x,sigma_HR,b_HR),1e-6,w)$value/w

############################################
# VREHT MODEL
############################################

  loglik_VR <- function(par){

    sigma <- par[1]
    a <- par[2]

    p <- interval_prob(function(x,p) g_VREHT(x,p[1],p[2]),c(sigma,a))

    -sum(counts*log(p))
  }

  fit_VR <- try(optim(c(25,5),loglik_VR,
                      method="L-BFGS-B",
                      lower=c(0.1,0.1),
                      upper=c(100,20),
                      hessian=TRUE),silent=TRUE)

  if(inherits(fit_VR,"try-error")) next

  sigma_VR <- fit_VR$par[1]
  a_VR <- fit_VR$par[2]

############################################
# AIC
############################################

  AIC_HN <- 2*fit_HN$value + 2
  AIC_HR <- 2*fit_HR$value + 4
  AIC_VR <- 2*fit_VR$value + 4

############################################
# STORE RESULTS
############################################

  valid_sim <- valid_sim + 1

  results <- rbind(results,data.frame(

    sim=valid_sim,

    pa_HN,pa_HR,pa_VR,

    sigma_HN,sigma_HR,b_HR,sigma_VR,a_VR,

    AIC_HN,AIC_HR,AIC_VR
  ))

  cat("Valid simulation:",valid_sim,"\n")
}

cat("Total attempts:",sim_attempt,"\n")

############################################
# MODEL SELECTION METRICS
############################################

AIC_mat <- cbind(results$AIC_HN,
                 results$AIC_HR,
                 results$AIC_VR)

deltaAIC <- t(apply(AIC_mat,1,function(x)x-min(x)))

PMS_HN <- mean(deltaAIC[,1] <=2)
PMS_HR <- mean(deltaAIC[,2] <=2)
PMS_VR <- mean(deltaAIC[,3] <=2)

best_model <- apply(AIC_mat,1,which.min)

SF_HN <- mean(best_model==1)
SF_HR <- mean(best_model==2)
SF_VR <- mean(best_model==3)

############################################
# SAVE OUTPUT
############################################

write.csv(results,"results_HN_HR_mix.csv",row.names=FALSE)

print("Simulation finished")