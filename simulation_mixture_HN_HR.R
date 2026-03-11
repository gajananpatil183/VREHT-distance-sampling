############################################################
# VREHT Detection Function Simulation Study
# True detection model: VREHT(sigma = 20, a = 5)
# Competing models: Half-Normal (HN), Hazard-Rate (HR), VREHT
# Data type: Grouped line-transect distance sampling
# Author: Gajanan Patil
############################################################

rm(list=ls())

############################################
## Simulation Settings
############################################

w <- 50                 # truncation distance
sigma_true <- 20         # true VREHT scale parameter
a_true <- 5              # true VREHT shape parameter
nsim <- 500              # number of simulation replicates
n <- 200                 # sample size
eps <- 1e-5              # step size for numerical derivatives

breaks <- seq(0,w,length.out=6)  # grouped distance intervals

############################################
# DETECTION FUNCTIONS
############################################

g_HN <- function(x,sigma){
  exp(-(x^2)/(2*sigma^2))
}

g_HR <- function(x,sigma,b){
  1-exp(-(x/sigma)^(-b))
}

g_VREHT <- function(x,sigma,a){
  1-(tanh(x/sigma))^a
}

############################################
# TRUE Pa (VREHT MODEL)
############################################

true_pa <- integrate(function(x) g_VREHT(x,sigma_true,a_true),1e-6,w)$value/w

############################################
# SIMULATE DISTANCES (VREHT TRUE MODEL)
############################################

simulate_VREHT <- function(n){

  y <- numeric()

  while(length(y)<n){

    x <- runif(n,0,w)
    u <- runif(n)

    accept <- u < g_VREHT(x,sigma_true,a_true)

    y <- c(y,x[accept])
  }

  y[1:n]
}

############################################
# INTERVAL PROBABILITIES
############################################

interval_prob <- function(gfun,param){

  k <- length(breaks)-1
  p <- numeric(k)

  denom <- try(integrate(function(x) gfun(x,param),
                         1e-6,w,
                         subdivisions=1000,
                         rel.tol=1e-6)$value,
               silent=TRUE)

  if(inherits(denom,"try-error")) return(rep(NA,k))

  for(i in 1:k){

    num <- try(integrate(function(x) gfun(x,param),
                         max(breaks[i],1e-6),
                         breaks[i+1],
                         subdivisions=1000,
                         rel.tol=1e-6)$value,
               silent=TRUE)

    if(inherits(num,"try-error")) return(rep(NA,k))

    p[i] <- num/denom
  }

  p
}

############################################
# STORAGE
############################################

results <- data.frame()

sim_attempt <- 0
valid_sim <- 0

############################################
# SIMULATION LOOP
############################################

while(valid_sim < nsim){

sim_attempt <- sim_attempt + 1

y <- simulate_VREHT(n)

grp <- cut(y,breaks,include.lowest=TRUE)
counts <- as.numeric(table(grp))

############################################
# HN MODEL
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

Pa_up <- integrate(function(x) g_HN(x,sigma_HN+eps),1e-6,w)$value/w
Pa_dn <- integrate(function(x) g_HN(x,sigma_HN-eps),1e-6,w)$value/w

grad <- (Pa_up-Pa_dn)/(2*eps)
se_pa_HN <- sqrt(grad^2 * cov_HN)

se_sigma_HN <- sqrt(diag(cov_HN))

############################################
# HR MODEL
############################################

loglik_HR <- function(par){

sigma <- par[1]
b <- par[2]

p <- interval_prob(function(x,p) g_HR(x,p[1],p[2]),c(sigma,b))

-sum(counts*log(p))
}

fit_HR <- try(optim(c(20,3),loglik_HR,hessian=TRUE),silent=TRUE)
if(inherits(fit_HR,"try-error")) next

sigma_HR <- fit_HR$par[1]
b_HR <- fit_HR$par[2]

cov_HR <- try(solve(fit_HR$hessian),silent=TRUE)
if(inherits(cov_HR,"try-error")) next

pa_HR <- integrate(function(x) g_HR(x,sigma_HR,b_HR),1e-6,w)$value/w

Pa_sigma_up <- integrate(function(x) g_HR(x,sigma_HR+eps,b_HR),1e-6,w)$value/w
Pa_sigma_dn <- integrate(function(x) g_HR(x,sigma_HR-eps,b_HR),1e-6,w)$value/w
dPa_dsigma <- (Pa_sigma_up-Pa_sigma_dn)/(2*eps)

Pa_b_up <- integrate(function(x) g_HR(x,sigma_HR,b_HR+eps),1e-6,w)$value/w
Pa_b_dn <- integrate(function(x) g_HR(x,sigma_HR,b_HR-eps),1e-6,w)$value/w
dPa_db <- (Pa_b_up-Pa_b_dn)/(2*eps)

grad <- c(dPa_dsigma,dPa_db)

se_pa_HR <- sqrt(t(grad) %*% cov_HR %*% grad)

se_sigma_HR <- sqrt(diag(cov_HR))[1]
se_b_HR <- sqrt(diag(cov_HR))[2]

############################################
# VREHT MODEL
############################################

loglik_VR <- function(par){

sigma <- par[1]
a <- par[2]

p <- interval_prob(function(x,p) g_VREHT(x,p[1],p[2]),c(sigma,a))

-sum(counts*log(p))
}

fit_VR <- try(optim(c(20,5),loglik_VR,hessian=TRUE),silent=TRUE)
if(inherits(fit_VR,"try-error")) next

sigma_VR <- fit_VR$par[1]
a_VR <- fit_VR$par[2]

cov_VR <- try(solve(fit_VR$hessian),silent=TRUE)
if(inherits(cov_VR,"try-error")) next

pa_VR <- integrate(function(x) g_VREHT(x,sigma_VR,a_VR),1e-6,w)$value/w

Pa_sigma_up <- integrate(function(x) g_VREHT(x,sigma_VR+eps,a_VR),1e-6,w)$value/w
Pa_sigma_dn <- integrate(function(x) g_VREHT(x,sigma_VR-eps,a_VR),1e-6,w)$value/w
dPa_dsigma <- (Pa_sigma_up-Pa_sigma_dn)/(2*eps)

Pa_a_up <- integrate(function(x) g_VREHT(x,sigma_VR,a_VR+eps),1e-6,w)$value/w
Pa_a_dn <- integrate(function(x) g_VREHT(x,sigma_VR,a_VR-eps),1e-6,w)$value/w
dPa_da <- (Pa_a_up-Pa_a_dn)/(2*eps)

grad <- c(dPa_dsigma,dPa_da)

se_pa_VR <- sqrt(t(grad) %*% cov_VR %*% grad)

se_sigma_VR <- sqrt(diag(cov_VR))[1]
se_a_VR <- sqrt(diag(cov_VR))[2]

############################################
# CHECK FAILED ITERATIONS
############################################

check_values <- c(
pa_HN,pa_HR,pa_VR,
se_pa_HN,se_pa_HR,se_pa_VR,
sigma_HN,sigma_HR,b_HR,sigma_VR,a_VR,
se_sigma_HN,se_sigma_HR,se_b_HR,se_sigma_VR,se_a_VR
)

if(any(!is.finite(check_values))) next

############################################
# AIC
############################################

AIC_HN <- 2*fit_HN$value + 2
AIC_HR <- 2*fit_HR$value + 4
AIC_VR <- 2*fit_VR$value + 4

############################################
# CHI-SQUARE TEST
############################################

expected_HN <- n*interval_prob(function(x,p) g_HN(x,p),sigma_HN)
expected_HR <- n*interval_prob(function(x,p) g_HR(x,p[1],p[2]),c(sigma_HR,b_HR))
expected_VR <- n*interval_prob(function(x,p) g_VREHT(x,p[1],p[2]),c(sigma_VR,a_VR))

chi_HN <- chisq.test(counts,p=expected_HN/sum(expected_HN))
chi_HR <- chisq.test(counts,p=expected_HR/sum(expected_HR))
chi_VR <- chisq.test(counts,p=expected_VR/sum(expected_VR))

############################################
# STORE VALID RESULT
############################################

valid_sim <- valid_sim + 1

results <- rbind(results,data.frame(

sim = valid_sim,

pa_HN,
pa_HR,
pa_VR,

se_pa_HN,
se_pa_HR,
se_pa_VR,

sigma_HN,
sigma_HR,
b_HR,
sigma_VR,
a_VR,

se_sigma_HN,
se_sigma_HR,
se_b_HR,
se_sigma_VR,
se_a_VR,

AIC_HN,
AIC_HR,
AIC_VR,

p_HN = chi_HN$p.value,
p_HR = chi_HR$p.value,
p_VR = chi_VR$p.value
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
# UNCONDITIONAL SUMMARY
############################################

table_uncond <- data.frame(

Model=c("HN","HR","VREHT"),

PMF=c(mean(results$p_HN>0.05),
      mean(results$p_HR>0.05),
      mean(results$p_VR>0.05)),

PMS=c(PMS_HN,PMS_HR,PMS_VR),

SF=c(SF_HN,SF_HR,SF_VR),

ARB=c(abs(mean(results$pa_HN)-true_pa)/true_pa,
      abs(mean(results$pa_HR)-true_pa)/true_pa,
      abs(mean(results$pa_VR)-true_pa)/true_pa),

NRMSE=c(sqrt(mean((results$pa_HN-true_pa)^2))/true_pa,
        sqrt(mean((results$pa_HR-true_pa)^2))/true_pa,
        sqrt(mean((results$pa_VR-true_pa)^2))/true_pa),

EmpSE_pa=c(sd(results$pa_HN),
           sd(results$pa_HR),
           sd(results$pa_VR)),

MeanSE_pa=c(mean(results$se_pa_HN),
            mean(results$se_pa_HR),
            mean(results$se_pa_VR)),

sigma_est=c(mean(results$sigma_HN),
            mean(results$sigma_HR),
            mean(results$sigma_VR)),

b_est=c(NA,
        mean(results$b_HR),
        NA),

a_est=c(NA,
        NA,
        mean(results$a_VR)),

sigma_SE=c(mean(results$se_sigma_HN),
           mean(results$se_sigma_HR),
           mean(results$se_sigma_VR)),

b_SE=c(NA,
       mean(results$se_b_HR),
       NA),

a_SE=c(NA,
       NA,
       mean(results$se_a_VR))
)

############################################
# CONDITIONAL SUMMARY
############################################

cond_HN <- subset(results,p_HN>0.05)
cond_HR <- subset(results,p_HR>0.05)
cond_VR <- subset(results,p_VR>0.05)

table_cond <- data.frame(

Model=c("HN","HR","VREHT"),

Chi_p=c(mean(cond_HN$p_HN),
        mean(cond_HR$p_HR),
        mean(cond_VR$p_VR)),

PMS=c(PMS_HN,PMS_HR,PMS_VR),

SF=c(SF_HN,SF_HR,SF_VR),

ARB=c(abs(mean(cond_HN$pa_HN)-true_pa)/true_pa,
      abs(mean(cond_HR$pa_HR)-true_pa)/true_pa,
      abs(mean(cond_VR$pa_VR)-true_pa)/true_pa),

NRMSE=c(sqrt(mean((cond_HN$pa_HN-true_pa)^2))/true_pa,
        sqrt(mean((cond_HR$pa_HR-true_pa)^2))/true_pa,
        sqrt(mean((cond_VR$pa_VR-true_pa)^2))/true_pa),

EmpSE_pa=c(sd(cond_HN$pa_HN),
           sd(cond_HR$pa_HR),
           sd(cond_VR$pa_VR)),

MeanSE_pa=c(mean(cond_HN$se_pa_HN),
            mean(cond_HR$se_pa_HR),
            mean(cond_VR$se_pa_VR)),

sigma_est=c(mean(cond_HN$sigma_HN),
            mean(cond_HR$sigma_HR),
            mean(cond_VR$sigma_VR)),

b_est=c(NA,
        mean(cond_HR$b_HR),
        NA),

a_est=c(NA,
        NA,
        mean(cond_VR$a_VR)),

sigma_SE=c(mean(cond_HN$se_sigma_HN),
           mean(cond_HR$se_sigma_HR),
           mean(cond_VR$se_sigma_VR)),

b_SE=c(NA,
       mean(cond_HR$se_b_HR),
       NA),

a_SE=c(NA,
       NA,
       mean(cond_VR$se_a_VR))
)

############################################
# PRINT RESULTS
############################################

write.csv(results,"simulation_results.csv",row.names=FALSE)
write.csv(table_uncond,"unconditional_summary.csv",row.names=FALSE)
write.csv(table_cond,"conditional_summary.csv",row.names=FALSE)

table_cond

