############################################################
# Simulation study for grouped line transect distance sampling
# Models: Half-Normal (HN), Hazard-rate (HR), VREHT
# Author: Gajanan Patil
# Description:
# This script simulates detection distances under mixture models,
# fits grouped detection functions, and computes performance metrics
# including Pa, SE(Pa), ARB, NRMSE, and AIC-based selection.
############################################################

rm(list=ls())
###Require Packages###
#install.packages(c("mrds","maxLik","Deriv","Distance"));install.packages("goftest");install.packages("DSsim");install.packages("nleqslv")
library(mrds);library(maxLik);library(Distance);library(Deriv);library(goftest);library(nleqslv)

################ SETTINGS ################
si=5     ###Put scale parametr of HN
a=25;b=3;  #####Put scale and shape paramter of HR

n=50 ## sample size 50,100,150,200 
w=50  ##Truncation disatnce
Phi1=0.5;Phi2=0.5  ## Coefficient of mixtute componant
nmin=10 ## Simulation sample
###########################################

####Intialize the names to store result
TNN1=0;TNN2=0;TNN3=0
AN1=0;AN2=0;AN3=0;SLA=0;pvHN=0;pvVHRT=0;pvHR=0;True_Pa=0;HN_pa=0;HR_pa=0;VHRT_pa=0;HN_AIC=0;HR_AIC=0;VHRT_AIC=0;HN_Chip=0;HR_Chip=0;VHRT_Chip=0
AN5=0;AN6=0;AN7=0;AN8=0;L_HN=0;U_HN=0;CL_HN=0;CP_HN=0;L_HR=0;U_HR=0;CL_HR=0;CP_HR=0;L_VHRT=0;U_VHRT=0;CL_VHRT=0;CP_VHRT=0
HN_se_sigma=0;HN_se_pa=0;HR_se_sigma=0;HR_se_b=0;HR_se_pa=0;VHRT_se_sigma=0;VHRT_se_a=0;VHRT_se_pa=0;ARB_Hn=0;ARB_Hr=0;ARB_VREHT=0
SLHn=0;SLHr=0;SLVreht=0;HN_si=0;HR_si=0;HR_b=0;VREHT_si=0;VREHT_a=0 

###########
E=1  
  while(E<nmin)
   {
     #### To generate random sample from True detection function model
     FF1=function(ZZ)  ######### to get the value of mu
         {
          ((Phi1)*(exp(-(ZZ^2)/(2*si^2))))+((Phi2)*(1-exp(-(ZZ/a)^(-b))))
          }
        	MUU=integrate(FF1,lower=0,upper=w)$value

	distance=0 ###### to store values of distance sample data
      v2=1
	while(v2<n+1) 
 		{
   		repeat{
   		     u=runif(1,0,1) 
  		     gp=function(x){
       	        fn1=function(q)
         	         {
                    ((Phi1)*(exp(-(q^2)/(2*si^2))))+((Phi2)*(1-exp(-(q/a)^(-b))))
          	         }
          	Tt=tryCatch(integrate(fn1,lower=0,upper=x),error=function(Tt){})
		    if(length(Tt)==0)
 		         {
                      dd=w+5
 		         }else{
   	 	          Tt=Tt
  	               }
       	    dd1=Tt$value-u*MUU
        	    return(dd1)
         	      }
            dd=tryCatch(nleqslv(.5,gp)$x,error=function(dd){})
        	break
         	}
              if(dd >w || length(dd)==0)
                {
                  next
                }else{
                  distance[v2]=dd                 
                } 
            v2=v2+1     
		}
       x=distance ##### This is the perpendicular distances
  #print(x)
##################To analyze this data##########

	k=5             # number of equal-width bins
	L=1000          # total transect length (replace with your survey length)

	edges = seq(0, w, length.out = k+1)
	if (exists("x") && !is.null(x)) {
  	counts = as.numeric(table(cut(x, breaks = edges, include.lowest = TRUE)))
	}
	n <- sum(counts)

      ### HN function
	g_hn = function(x, par) {
	sigma <- par[1]
  	exp(- x^2 / (2 * sigma^2))
	}

	##HR function
	g_hr = function(x, par) {
  	sigma <- par[1]; b <- par[2]
  	1 - exp(- (x / sigma)^(-b))
	}

	##VREHT function
	g_vreht <- function(x, par) {
  	sigma <- par[1]; a <- par[2]
  	1 - (tanh(x / sigma))^a
	}

	# helper: numeric integral of g from lo to hi
	int_g = function(gfun, par, lo, hi) {
  	f = function(xx) gfun(xx, par)
  	out = try(integrate(f, lower = lo, upper = hi, rel.tol = 1e-8, subdivisions = 1000), silent = TRUE)
  	if (inherits(out, "try-error")) {
    	warning("integration failed; returning NA")
    	return(NA_real_)
  	}
  	out$value
	}

	bin_probs_from_g = function(gfun, par, edges, w) {
  	nums = numeric(length(edges)-1)
  	for (i in seq_len(length(nums))) {
    	nums[i] = int_g(gfun, par, edges[i], edges[i+1])
  	}
  	denom = sum(nums)
  	if (denom <= 0 || any(is.na(nums))) {
    	return(rep(0, length(nums)))
  	}
  	nums / denom
	}

	nll_multinomial = function(logpar, gfun, counts, edges, w) {
  		par = exp(logpar)  # positivity
  		p = bin_probs_from_g(gfun, par, edges, w)
  		if (any(p <= 0) || any(is.na(p))) return(1e12)   # penalize invalid params
  		- sum(counts * log(p))
		}

	fit_grouped_model = function(gfun, init_par, counts, edges, w) {
  		init_logpar = log(init_par)
    		## Maximize Likelihood this using Otim function
  		opt = optim(par = init_logpar,fn = nll_multinomial,gfun = gfun,counts = counts,edges = edges,w = w, hessian = TRUE, method = "BFGS",control = list(maxit = 2000, reltol = 1e-8))
  		par_hat = exp(opt$par)
  		mu_hat = int_g(gfun, par_hat, 0, w)
  		pa_hat = mu_hat / w  ####Estiamte value of Pa

  		D_hat = if (!is.na(mu_hat) && mu_hat > 0) n / (2 * L * mu_hat) else NA_real_
  		fitted_p = bin_probs_from_g(gfun, par_hat, edges, w)
  		list(par_hat = par_hat, loglik = -opt$value, converged = (opt$convergence==0),
       	mu_hat = mu_hat, pa_hat = pa_hat, D_hat = D_hat, fitted_p = fitted_p, opt = opt)
		}

	init_hn=c(sigma = 20)       # HN: sigma
	init_hr=c(sigma = w/4, b = 2) # HR: sigma, b
	init_vre=c(sigma =20, a = 2) # VREHT: sigma, a

	fit_hn=fit_grouped_model(g_hn,  init_hn,  counts, edges, w)
	fit_hr=fit_grouped_model(g_hr,  init_hr,  counts, edges, w)
	fit_vre=fit_grouped_model(g_vreht, init_vre, counts, edges, w)
	
######Since the simulated data arise from diverse detection shapes, numerical integration and parameter estimation occasionally failed to converge or produced 
       #####unrealistic parameter values. To ensure numerical stability and reliability of the results, such iterations were identified based on predefined parameter 
         ###thresholds and subsequently excluded from the analysis.
	if((fit_vre$par_hat[2] > 1e3) ||(fit_vre$par_hat[1] > 1e3)|| (fit_hr$par_hat[1]<(1e-4))||(fit_vre$par_hat[2] < 1e-3) || (fit_hr$par_hat[1]>(1e4)))
	{
	E=E
	}else{
	aic_hn5  <- -2 * fit_hn$loglik + 2 * 1 ###AIC of HN
	aic_hr5  <- -2 * fit_hr$loglik + 2 * 2 ### AIC of HR
	aic_vre5 <- -2 * fit_vre$loglik + 2 * 2 ### AIC of VREHT

 ############# Density and SD estiamte of HN  ##########

	mu_Hn=fit_hn$mu_hat
	pa_Hn=mu_Hn/w            ###Pa estimate using HN
	AIC_Hn=aic_hn5

	par_hat = fit_hn$par_hat
	opt = fit_hn$opt

	# --- parameter SE ---
	vcov_log = solve(opt$hessian)
      vcov_log1=tryCatch(vcov_log = solve(opt$hessian), error=function(e) NULL)
       if(length(vcov_log1)==0)
            {
		vcov_log=1
		}else{
		vcov_log=vcov_log1
           }
	J1t = diag(par_hat, nrow=1, ncol=1)
	vcov_theta = J1t %*% vcov_log %*% J1t
	se_par = sqrt(diag(vcov_theta))

	# --- pa function ---
	pa_fun <- function(par){
 	integrate(function(x) exp(-(x^2)/(2*par[1]^2)),0,w)$value / w
	}

	# --- gradient ---
	grad_pa_fun <- function(par, eps=1e-6){
 	ghn <- numeric(length(par))
 	for(j in 1:length(par)){
  	p1=par; p2=par
  	p1[j]=p1[j]+eps
  	p2[j]=p2[j]-eps
  	ghn[j]=(pa_fun(p1)-pa_fun(p2))/(2*eps)
 	}
 	ghn
	}
	grad = grad_pa_fun(par_hat)

	# --- SE(pa) ---
	var_pa = t(grad) %*% vcov_theta %*% grad #####varaicne pf Pa using HN


############# Density and SD estiamte of HR  ##########

	mu_Hr=fit_hr$mu_hat
	pa_Hr=mu_Hr/w   #####Pa using HR
	AIC_Hr=aic_hr5

	# ===== HR SE of paramter and Pa=====
	par_hathr = fit_hr$par_hat
	opt_hr = fit_hr$opt

	# --- parameter SE ---
	k_hr = length(par_hathr)
      vcov_loghr1=tryCatch(solve(opt_hr$hessian), error=function(e) NULL)
       if(length(vcov_loghr1)==0)
            {
		vcov_loghr=diag(par_hathr, nrow = k_hr, ncol = k_hr)
		}else{
		vcov_loghr=vcov_loghr1
           }
	Jhr = diag(par_hathr, nrow = k_hr, ncol = k_hr)
	vcov_thetahr = t(Jhr) %*% vcov_loghr %*% Jhr
	se_parhr = sqrt(diag(vcov_thetahr))

	# --- pa function ---
	pa_fun_hr = function(par){
 	integrate(function(x) 1-exp(-(x/par[1])^(-par[2])),0,w)$value / w
	}

	# --- gradient ---
	grad_pa_hr = function(par, eps=1e-6){
 	ghr<- numeric(length(par))
 	for(j in 1:length(par)){
  		p1=par; p2=par
  		p1[j]=p1[j]+eps
  		p2[j]=p2[j]-eps
  		ghr[j]=(pa_fun_hr(p1)-pa_fun_hr(p2))/(2*eps)
 		}
 	ghr
	}
	gradhr = grad_pa_hr(par_hathr)

	# --- SE(pa) ---
	var_pahr = t(gradhr) %*% vcov_thetahr %*% gradhr   ###Varaicne of Pa using HR


############### Density and SD estiamte of VREHT###############

	mu_TN=fit_vre$mu_hat
	pa_TN=mu_TN/w     ### Pa estimate using VREHT
	AIC_TN=aic_vre5

	# ===== VREHT =====
	par_hatvre = fit_vre$par_hat
	opt_vre = fit_vre$opt

	# --- parameter SE ---
	k_vre=length(par_hatvre)
	#vcov_logvre = solve(opt_vre$hessian)
      vcov_logvre1=tryCatch(solve(opt_vre$hessian), error=function(e) NULL)
       if(length(vcov_logvre1)==0)
            {
		vcov_logvre=diag(par_hatvre, nrow = k_vre, ncol = k_vre)
		}else{
		vcov_logvre=vcov_logvre1
           }
	Jvre = diag(par_hatvre, nrow = k_vre, ncol = k_vre)
	vcov_thetavre = t(Jvre) %*% vcov_logvre %*% Jvre
	se_parvre = sqrt(diag(vcov_thetavre))

	# --- pa function ---
	pa_fun_vre = function(par){
 	integrate(function(x) 1-(tanh(x/par[1]))^par[2],0,w)$value / w
	}

	# --- gradient ---
	grad_pa_vre = function(par, eps=1e-6){
 	gvr <- numeric(length(par))
 	for(j in 1:length(par)){
  		p1=par; p2=par
  		p1[j]=p1[j]+eps
  		p2[j]=p2[j]-eps
  		gvr[j]=(pa_fun_vre(p1)-pa_fun_vre(p2))/(2*eps)
 		}
 		gvr
		}
	gradvre = grad_pa_vre(par_hatvre)

	# --- SE(pa) ---
	var_pavre = t(gradvre) %*% vcov_thetavre %*% gradvre ## Varaicne of Pa using VREHT

# -----------------------------------------------------------
# Pearson chi-square GOF test for grouped multinomial models

gof_chisq=function(fit, counts, n, k, q) {
  expected=n * fit$fitted_p
  chisq=sum((counts - expected)^2 / expected)
  df=k - 1 - q
  pval=1 - pchisq(chisq, df)
  list(ChiSq = chisq, df = df, p.value = pval,observed = counts, expected = round(expected,2))
}

# Run GOF for each model
gof_hn=gof_chisq(fit_hn,  counts, n, k, q = 1)   
gof_hr=gof_chisq(fit_hr,  counts, n, k, q = 2)  
gof_vre=gof_chisq(fit_vre, counts, n, k, q = 2) 

###### chisq test in line transect sampling VREHT##########
		CVM_PvalhnL=gof_hn$p.value  ###Chi square P value of HN
	      CVM_Pvalhr=gof_hr$p.value   ###Chi square P value of HR



      CVM_PvalTn=gof_vre$p.value       ###Chi square P value of VREHT


  if( (fit_vre$opt$counts[2]>4)||(length(vcov_log1)=0) ||(length(vcov_loghr1)==0) || (length(vcov_logvre1)==0) )
     {
      E=E
     }else{

      
#############################Store all values 
      True_Pa=MUU/w     ## True value of Pa
      HN_pa[E]=pa_Hn    ## Pa estiamte using HN
      HR_pa[E]=pa_Hr	## Pa estiamte using HN
      VHRT_pa[E]=pa_TN  ## Pa estiamte using VREHT

      ### Absolute relative bias of HN, HR and VREHT
	ARB_Hn[E]=abs(pa_Hn-True_Pa)
	ARB_Hr[E]=abs(pa_Hr-True_Pa)
	ARB_VREHT[E]=abs(pa_TN-True_Pa)

	## Parameter estimate of HN, HR and VREHT 
	HN_si[E] =fit_hn$par_hat[1]
	HR_si[E] = fit_hr$par_hat[1]
	HR_b[E]= fit_hr$par_hat[2]
	VREHT_si[E] =fit_vre$par_hat[1]
	VREHT_a[E]=fit_vre$par_hat[2]

	##SE of parameter estiamte of HN, HR and VREHT
	HN_se_sigma[E] = se_par[1]
	HR_se_sigma[E] = se_parhr[1]
	HR_se_b[E]= se_parhr[2]
	VHRT_se_sigma[E] = se_parvre[1]
	VHRT_se_a[E]= se_parvre[2]

	###SE of Pa estiamte of HN< HR and VREHT
	HN_se_pa[E] = sqrt(var_pa)
	HR_se_pa[E] = sqrt(var_pahr)
	VHRT_se_pa[E] = sqrt(var_pavre)

	##AIC of HN,HR and VREHT
      HN_AIC[E]=AIC_Hn
      HR_AIC[E]=AIC_Hr
      VHRT_AIC[E]=AIC_TN
 
	### Chi square goodness of fit test p value of HN, HR and VREHT 
      HN_Chip[E]=CVM_PvalhnL
      HR_Chip[E]=CVM_Pvalhr
      VHRT_Chip[E]=CVM_PvalTn

	###Counting number of model fit     
      pvHN[E]=ifelse(CVM_PvalhnL<=0.05,0,1)
      pvHR[E]=ifelse(CVM_Pvalhr<=0.05,0,1)
      pvVHRT[E]=ifelse(CVM_PvalTn<=0.05,0,1)

     ########## to find delta AIC#####
     an1=AIC_Hn-min(AIC_Hn,AIC_Hr,AIC_TN)
     an2=AIC_Hr-min(AIC_Hn,AIC_Hr,AIC_TN)
     an3=AIC_TN-min(AIC_Hn,AIC_Hr,AIC_TN)
     
    AN6[E]=ifelse((an1<=2) && (CVM_PvalhnL > 0.05) ,1,0);AN7[E]=ifelse((an2<=2) &&(CVM_Pvalhr > 0.05),1,0);AN8[E]=ifelse(an3<=2 && (CVM_PvalTn> 0.05),1,0)     

	SLHn[E]=ifelse(AIC_Hn==min(AIC_Hn,AIC_Hr,AIC_TN),1,0)
	SLHr[E]=ifelse(AIC_Hr==min(AIC_Hn,AIC_Hr,AIC_TN),1,0)
	SLVreht[E]=ifelse(AIC_TN==min(AIC_Hn,AIC_Hr,AIC_TN),1,0)

  #### Model selectiob metrix using delta AIC
   p1hn=CVM_PvalhnL;p1hr=CVM_Pvalhr;p1tn=CVM_PvalTn;  

  if(p1hn>0.05)
     {
       if(p1hr >0.05 && p1tn >0.05 && (AIC_Hn <= min(AIC_Hr,AIC_TN))){SS=1}
       if(p1hr >0.05 && p1tn >0.05 && (AIC_Hr <= min(AIC_Hn,AIC_TN))){SS=2}
       if(p1hr >0.05 && p1tn >0.05 && (AIC_TN <= min(AIC_Hr,AIC_Hn))){SS=3}
       
       if(p1hr >0.05 && p1tn <= 0.05 && (AIC_Hn <=(AIC_Hr))){SS=1}
       if(p1hr >0.05 && p1tn <= 0.05 && (AIC_Hr <=(AIC_Hn))){SS=2}

       if(p1hr <= 0.05 && p1tn > 0.05 && (AIC_Hn <=(AIC_TN))){SS=1}
       if(p1hr <= 0.05 && p1tn > 0.05 && (AIC_TN <=(AIC_Hn))){SS=3}
       if(p1hr <= 0.05 && p1tn <=0.05) {SS=1}

      }else{
      if(p1hr > 0.05 && p1tn >  0.05 && (AIC_Hr <= AIC_TN)){SS=2}
      if(p1hr > 0.05 && p1tn >  0.05 && (AIC_TN <= AIC_Hr)){SS=3}
      if(p1hr > 0.05 && p1tn <= 0.05) {SS=2}
      if(p1hr <= 0.05 && p1tn > 0.05) {SS=3}
      if(p1hr <= 0.05 && p1tn <= 0.05 ){SS=0}
     }
    SLA[E]=SS
	
        E=E+1
  	}
	}
	  next	     	   
   }

mg=length(SLA)
#### Store all data
ut=data.frame(SLA,True_Pa,HN_pa,HR_pa,VHRT_pa,HN_AIC,HR_AIC,VHRT_AIC,HN_Chip,HR_Chip,VHRT_Chip,n,mg,ARB_Hn,ARB_Hr,ARB_VREHT,HN_se_pa,HR_se_pa,VHRT_se_pa,HN_si,HR_si,HR_b,VREHT_si,VREHT_a,HN_se_sigma,HR_se_sigma,HR_se_b,VHRT_se_sigma,VHRT_se_a,SLHn,SLHr,SLVreht,AN6,AN7,AN8)         
	

ut2=ut[order(SLA),]
gg5=colMeans(ut2)
MShn=length(which(SLA==1))
MShr=length(which(SLA==2))
MSvhrt=length(which(SLA==3))

#### To store uncoditional result
ut88=data.frame(n,mg,"TRUE_Pa"=mean(True_Pa),"SL_HN"=(sum(SLHn)*100)/mg,"SLHR"=(sum(SLHr)*100)/mg,"SLVREHT"=(sum(SLVreht)*100)/mg,"pahn"=mean(HN_pa),"pahr"=mean(HR_pa),"paVHRT"=mean(VHRT_pa),"ARB_HN"=mean(ARB_Hn)/True_Pa,"ARB_HR"=mean(ARB_Hr)/True_Pa,"ARB_VREHT"=mean(ARB_VREHT)/True_Pa,"MES_HN"=((((mean(HN_pa)-mean(True_Pa))^2)+var(HN_pa)))/True_Pa,"MSE_HR"=((((mean(HR_pa)-mean(True_Pa))^2)+var(HR_pa)))/True_Pa,"MSE_VHRT"=((((mean(VHRT_pa)-mean(True_Pa))^2)+var(VHRT_pa)))/True_Pa,"EpSEHN"=sqrt(var(HN_pa)),"EpSEHR"=sqrt(var(HR_pa)),"EpSEVREHT"=sqrt(var(VHRT_pa)),"SEpa HN"=mean(HN_se_pa),"SEpa HR"=mean(HR_se_pa),"SEpa"=mean(VHRT_se_pa),"HN_si"=mean(HN_si),"HR_si"=mean(HR_si),"HR_b"=mean(HR_b),"VREHT_si"=mean(VREHT_si),"VREHT_a"=mean(VREHT_a),"HN_se_sigma"=mean(HN_se_sigma),"HR_se_sigma"=mean(HR_se_sigma),"HR_se_b"=mean(HR_se_b),"VHRT_se_sigma"=mean(VHRT_se_sigma),"VHRT_se_a"=mean(VHRT_se_a),"HN fit"=(sum(pvHN)*100/mg),"HR fit"=(sum(pvHR)*100/mg),"VHRT fit"=(sum(pvVHRT)*100/mg))                 

##### To store conditonal (fitted model) Result

Report_HN=data.frame(True_Pa,HN_pa)
Report_HR=data.frame(True_Pa,HR_pa)
Report_VRHT=data.frame(True_Pa,VHRT_pa)
DaicHN=mean(AN6)*100;DaicHR=mean(AN7)*100;DaicVRHT=mean(AN8)*100;

ut71=ut[which(HN_Chip>=0.05),c(1,9,3,14,17,20,25,30)]
mse_hn1=var(ut71$HN_pa)+(mean(ut71$HN_pa-True_Pa))^2
HN_rep1=c(colMeans(ut71),"NARB"=mean(ut71$ARB_Hn/True_Pa),"NRMSE_Pa"=mse_hn1/True_Pa,"EpSEhn"=sqrt(var(ut71$HN_pa)),"DaicHN"=DaicHN)
HN_rep <- c(length(SLA),True_Pa,HN_rep1[1:6],emty1 = 0,HN_rep1[7],emty2 = 0,HN_rep1[8:length(HN_rep1)])

ut72=ut[which(HR_Chip>=0.05),c(1,10,4,15,18,21,22,26,27,31)]
mse_hr1=(var(ut72$HR_pa)+(mean(ut72$HR_pa-True_Pa))^2)
HR_rep=c(length(SLA),True_Pa,colMeans(ut72),"NARB"=mean(ut72$ARB_Hr/True_Pa),"NRMSE_Pa"=mse_hr1/True_Pa,"EpSEhr"=sqrt(var(ut72$HR_pa)),"DaicHR"=DaicHR)

ut73=ut[which(VHRT_Chip>=0.05),c(1,11,5,16,19,23,24,28,29,32)]
mse_vr1=(var(ut73$VHRT_pa)+(mean(ut73$VHRT_pa-True_Pa))^2)
VREHT_rep=c(length(SLA),True_Pa,colMeans(ut73),"NARB"=mean(ut73$ARB_VREHT/True_Pa),"NRMSE_Pa"=mse_vr1/True_Pa,"EpSEVREHT"=sqrt(var(ut73$VHRT_pa)),"DaicVREHT"=DaicVRHT)

ut74=c("total sample","True pa","SLA","Chi p","Pa_hat","ave ARB","SE_pa","sc par","sp par","SE sc","SE sp","SF","NARB","NRMSE","Emp SE pa","D AIC")
length(VREHT_rep)
ut88
ut75=data.frame(ut74,HN_rep,HR_rep,VREHT_rep)
ut75
write.csv(ut88,"simulation_unconditional_results.csv",row.names=FALSE)
write.csv(ut75,"simulation_conditional_results.csv",row.names=FALSE)
