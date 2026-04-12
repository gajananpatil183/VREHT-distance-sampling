# ============================================================
# MCDS Simulation (Point-transect): Covariate-based Detection Functions
# Author: Gajanan Patil and Shashibhushan Mahadik
# Description:
# This script implements a Monte Carlo simulation study under the 
# Multiple Covariate Distance Sampling (MCDS) framework for point transect data.
#
# Radial distance data are generated using a mixture detection function
# with both continuous and categorical covariates. The scale parameter of the
# detection function is modeled as a function of covariates using a log-link.
#
# Three detection models are fitted and compared:
# Half-Normal (HN), Hazard-Rate (HR), and the proposed VREHT model.
#
# Model performance is evaluated using Akaike Information Criterion (AIC),
# Cramér–von Mises goodness-of-fit test, detection probability (Pa),
# density estimates (D), bias, mean squared error (MSE), and confidence interval length (CIL).
#
# The results are summarized to assess the effect of covariates on detection
# probability and to compare the efficiency of competing models.
# ============================================================

# ---------------------------
# Initialization
# ---------------------------
rm(list = ls())

# ---------------------------
# Required Libraries
# ---------------------------
library(maxLik)
library(DSsim)
library(mrds)
library(Deriv)
library(Rdistance)
library(fastDummies)
library(stringi)
library(data.table)
library(rlang)
library(vctrs)

# ---------------------------
# Mixture Model Coefficients
# ---------------------------
Phi1 <- 0.5
Phi2 <- 0.5

# ---------------------------
# Survey Design (Point Transect)
# ---------------------------
w  <- 50     # Truncation distance (radius)
Ns <- 100    # Initial sample size

# Number of points (derived)
K1 <- round((Ns / (w * w * pi)))
K  <- ifelse(K1 == 0, 1, K1)

# ---------------------------
# Covariate Model Parameters
# ---------------------------

## Half-Normal (HN)
Bo0 <- 1.5
Bo1 <- 0.5
B21 <- 0.4
B22 <- 0.6

## Hazard-Rate (HR)
BB0  <- 2.5
BB1  <- 1
BB21 <- 0.8
BB22 <- 1

AA0 <- 1   # Shape parameter

# ---------------------------
# Storage Initialization
# ---------------------------

MM <- 0

# Density estimates
D_Hn1    <- 0
D_HR1    <- 0
D_vreht1 <- 0

# Density (alternative form)
D_Hn19    <- 0
D_HR19    <- 0
D_vreht19 <- 0

# Detection probability (Pa)
piHn3 <- 0
piHr3 <- 0
piVr3 <- 0

# Goodness-of-fit (CvM p-values)
pvalHn   <- 0
pvalHr   <- 0
pvalVRHT <- 0

# Model fit indicators
pvHN   <- 0
pvHR   <- 0
pvVHRT <- 0

# Model selection indicators
AN6 <- 0
AN7 <- 0
AN8 <- 0

# AIC values
AIC_Hn1 <- 0
AIC_Hr1 <- 0
AIC_TN1 <- 0

nn=501 #simulation size

J <- 1
while(J<nn)
{
op8=runif(Ns,0,1)
op=w*sqrt(op8)
ind1=which(op<=w)
xx1=op[ind1]
z1x=round(rchisq(Ns,3))

namelist<-c("A", "B", "C")
g50=sample(namelist,Ns,replace=TRUE) 

g51=data.frame(z1x,g50)
g52=dummy_cols(g51, select_columns = c("g50"),remove_first_dummy = TRUE)


z21=g52$g50_B;z22=g52$g50_C
mx=length(xx1);
Ara=300


  muux=0
   for( j in 1:mx)
  {
     f3x=function(z)
                {					
			Phi1*z*exp(-(z^2)/(2*(exp(Bo0+Bo1*z1x[j]+B21*z21[j]+B22*z22[j])^2)))+Phi2*z*(1-exp(-(z/(exp(BB0+BB1*z1x[j]+BB21*z21[j]+BB22*z22[j])))^(-AA0)))	
                  }
               ttx=integrate(f3x,lower=0,upper=w)
               muux[j]=ttx$value
  }  
  h0TnPM=1/muux
  pihnx=(2/w^2)*(1/h0TnPM)          
  True_Pa=mean(pihnx)

u1=runif(Ns,0,1);u2=ifelse(u1<=pihnx,1,0)
d=data.frame(xx1,z1x,z21,z22,g50,pihnx,u1,u2)
head(d)
d2=d[which(u2==1),]
head(d2)
m=length(d2$xx1);z1=d2$z1x;x=d2$xx1;Z21=d2$z21;Z22=d2$z22

lm=mean(z1)

###################### VREHT #############################
 f=0;mu=0;
 qq=function(par)
     {
	 b0=par[1]
       b1=par[2]
       a=par[3];b21=par[4];b22=par[5]
       for(k in 1:m)
         {
           f2o=function(z)
               {
		zz=z*(1-(tanh(z/(exp(b0+b1*z1[k]+b21*Z21[k]+b22*Z22[k]))))^exp(a))		
 		#zz1=ifelse(zz==NaN,1e-300,zz)
             return(zz)
               }
            tt=tryCatch(integrate(f2o,lower=0,upper=w,stop.on.error=FALSE),error=function(tthr1){})
            if(length(tt)==0)
             {
                       mu[k]=0.5
                       M49=0
               }else{               
                  mu[k]=tt$value
                 }              
 		f[k]=(x[k])*(1-(tanh(x[k]/(exp(b0+b1*z1[k]+b21*Z21[k]+b22*Z22[k]))))^exp(a))/mu[k]
            ifelse(f[k]<1e-300,1e-300,f[k])        
	    }
	  ww=prod(f) 
        ww1=ifelse(ww<1e-300,1e-300,ww)
	  f3=log(ww1)   
        return(f3)
      }
res12 <- maxBFGS(qq, start=c(5,(1/lm),1,.5,.5))
res=tryCatch(res12,error=function(res12){})
summary(res)


########### HN ###############################################
muhn=0;fhn=0
qqhn=function(par)
     {
	 b0=par[1]
       b1=par[2]
       bh21=par[3];bh22=par[4]
       for(k in 1:m)
         {
           f2hn=function(z)
               {
            zzhn=z*exp(-(z^2)/(2*(exp(b0+b1*z1[k]+bh21*Z21[k]+bh22*Z22[k])^2)))			
 		#zz1=ifelse(zz==NaN,1e-300,zz)
             return(zzhn)
               }
            tthn1=integrate(f2hn,lower=0,upper=w,stop.on.error=FALSE)
            muhn[k]=tthn1$value
 		fhn[k]=x[k]*(exp(-(x[k]^2)/(2*(exp(b0+b1*z1[k]+bh21*Z21[k]+bh22*Z22[k])^2))))/muhn[k]
            ifelse(fhn[k]<1e-300,1e-300,fhn[k])        
	    }
	  wwhn=prod(fhn) 
        ww1hn=ifelse(wwhn<1e-300,1e-300,wwhn)
	  f3hn=log(ww1hn)   
        return(f3hn)
      }
resHn12<- maxBFGS(qqhn, start=c(5,(1/lm),0.5,0.5))
resHn=tryCatch(resHn12,error=function(resHn12){})
summary(resHn)

############ HR ###############
fhr=0;muhr=0;
 qqhr=function(par)
     {
	 b0hr=par[1]
       b1hr=par[2]
       ahr=par[3];bhr21=par[4];bhr22=par[5]
       for(k in 1:m)
         {
           f2hr=function(z)
               {
		zzhr=z*(1-exp(-(z/(exp(b0hr+b1hr*z1[k]+bhr21*Z21[k]+bhr22*Z22[k])))^(-exp(ahr))))		
             return(zzhr)
               }
            ttHr=integrate(f2hr,lower=0,upper=w,stop.on.error=FALSE)
            muhr[k]=ttHr$value
 		fhr[k]=x[k]*(1-exp(-(x[k]/(exp(b0hr+b1hr*z1[k]+bhr21*Z21[k]+bhr22*Z22[k])))^(-exp(ahr))))/muhr[k]
            ifelse(fhr[k]<1e-300,1e-300,fhr[k])        
	    }
	  wwhr=prod(fhr) 
        ww1hr=ifelse(wwhr<1e-300,1e-300,wwhr)
	  f3hr=log(ww1hr)   
        return(f3hr)
      }
resHr12 <- maxBFGS(qqhr, start=c(3,(1/lm),1,.2,.2))
resHr=tryCatch(resHr12,error=function(resHr12){})
summary(resHr)
###############################################################################
 if((M49=0  || length(res)==0) || (length(resHn)==0)|| (length(resHr)==0)|| is.na(res$gradient[1]) || is.na(res$gradient[2]) || is.na(res$gradient[3]) || is.na(resHr$gradient[1]) || is.na(resHr$gradient[2]) || is.na(resHr$gradient[3]) || (resHr$gradient[1]==0 && resHr$gradient[2]==0)||(-0.01 >= res$gradient[1] || res$gradient[1] >= 0.01 ) || (-0.01 >= res$gradient[2] || res$gradient[2] >= 0.01 ) || (-0.01 >= res$gradient[3] || res$gradient[3] >= 0.01 ) || (-0.01 >= resHr$gradient[1] || resHr$gradient[1] >= 0.01 ) || (-0.01 >= resHr$gradient[2] || resHr$gradient[2] >= 0.01 ) || (-0.01 >= resHr$gradient[3] || resHr$gradient[3] >= 0.01 ))
   {
     J=J
       }else{
##VREHT##
#print(summary(res))

  B0=res$estimate[1]
  B1=res$estimate[2]
  AA=exp(res$estimate[3])
  BB12=res$estimate[4];BB22=res$estimate[5]
  muu=0
   	for( j in 1:m)
  	{
     	f3=function(z)
               {	
			z*(1-(tanh(z/(exp(B0+B1*z1[j]+BB21*Z21[j]+BB22*Z22[j]))))^AA)						  
                  }
               tt=integrate(f3,lower=0,upper=w)
               muu[j]=tt$value
  	}  

      h0_hat61=1/muu
      piTnL=mean(2/(w*w*h0_hat61))             
	D_vreht=(1/mean(piTnL))
      D_vreht18=(sum(h0_hat61))/(2*K*pi)

	piTnL12=mean(piTnL)

##########HN

  B0Hn=resHn$estimate[1]
  B1Hn=resHn$estimate[2]
  B21Hn=resHn$estimate[3]
  B22Hn=resHn$estimate[4]
  MUhn=0
   for( j in 1:m)
  {
     f3hn=function(z)
                {					
			z*exp(-(z^2)/(2*(exp(B0Hn+B1Hn*z1[j]+B21Hn*Z21[j]+B22Hn*Z22[j])^2)))	  
                  }
               ttHn=integrate(f3hn,lower=0,upper=w)
               MUhn[j]=ttHn$value
  }   
    

      h0_hat6=1/MUhn
      piHn=mean(2/(w*w*h0_hat6))       
	D_Hn=1/(mean(piHn))
      D_Hn18=(sum(h0_hat6))/(2*K*pi)
	piHn12=mean(piHn)

	
#########HR

  	B0Hr=resHr$estimate[1]
  	B1Hr=resHr$estimate[2]
  	AAHr=exp(resHr$estimate[3])
  	B21Hr=resHr$estimate[4]
  	B22Hr=resHr$estimate[5]

  muuHr=0
   for( j in 1:m)
  {
     f3HR=function(z)
               {	
			z*(1-exp(-(z/(exp(B0Hr+B1Hr*z1[j]+B21Hr*Z21[j]+B22Hr*Z22[j])))^(-AAHr)))						  
                  }
               ttHR=integrate(f3HR,lower=0,upper=w)
               muuHr[j]=ttHR$value
  }            
 	h0_hat62=1/muuHr
      piHr=mean(2/(w*w*h0_hat62))             
	D_HR=(1/(mean(piHr)))
      D_HR18=(sum(h0_hat62))/(2*K*pi)
	piHr12=mean(piHr)

###################
####################################################################################################
######HN CvM######

MUu5=0
  for( e in 1:m)
  {
     fx3=function(z)
               {					
		     z*exp(-(z^2)/(2*(exp(B0Hn+B1Hn*z1[e]+B21Hn*Z21[e]+B22Hn*Z22[e])^2)))	  		
                  }
               tx5=integrate(fx3,lower=0,upper=w)
               MUu5[e]=tx5$value
  }            
 #MUu
 tx25=0
 for(h in 1:m)
 {
   fk1=function(z25)
                 {
		     z25*exp(-(z25^2)/(2*(exp(B0Hn+B1Hn*z1[h]+B21Hn*Z21[h]+B22Hn*Z22[h])^2)))	  		
                  }
   tx25[h]=integrate(fk1,lower=0,upper=x[h])$value
  }
  cdf_hno5=tx25/MUu5

#######CVM
tada1=data.frame(cdf_hno5,x,z1,Z21,Z22)
tada=tada1[order(tada1$cdf_hno5), ]
library('goftest')
F_sto=0
ttt=tada$x
Zz11=tada$z1;
zz21=tada$Z21;
zz22=tada$Z22;
head(tada1)
n=m
 MUu22=0
  for( e in 1:m)
  {
     fx3=function(z)
               {					
		     z*exp(-(z^2)/(2*(exp(B0Hn+B1Hn*Zz11[e]+B21Hn*zz21[e]+B22Hn*zz22[e])^2)))	  		
                  }
               tx3=integrate(fx3,lower=0,upper=w)
               MUu22[e]=tx3$value
  }            
 tx2=0
 for(h in 1:m)
 {
   fk13=function(z22)
                 {
		     z22*exp(-(z22^2)/(2*(exp(B0Hn+B1Hn*Zz11[h]+B21Hn*zz21[h]+B22Hn*zz22[h])^2)))	  		
                  }
   tx2[h]=integrate(fk13,lower=0,upper=ttt[h])$value
  }
  cdf_hno3=tx2/MUu22
model_Hn = 1/(12*m) + sum((cdf_hno3 -((1:m)-.5)/m)^2)
CVM_Pval_Hn=1-pCvM(model_Hn, n = m)

################
######HR CvM###########################################################
#######################################################################
MUu7=0
  for( e in 1:m)
  {
     fx37=function(z)
               {					
			z*(1-exp(-(z/(exp(B0Hr+B1Hr*z1[e]+B21Hr*Z21[e]+B22Hr*Z22[e])))^(-AAHr)))						  
                  }
               tx7=integrate(fx37,lower=0,upper=w)
               MUu7[e]=tx7$value
  }            
 tx27=0
 for(h in 1:m)
 {
   fk17=function(z27)
                 {
			z27*(1-exp(-(z27/(exp(B0Hr+B1Hr*z1[h]+B21Hr*Z21[h]+B22Hr*Z22[h])))^(-AAHr)))						  
                  }
   tx27[h]=integrate(fk17,lower=0,upper=x[h])$value
  }
  cdf_hr7=tx27/MUu7

#######CVM
tada7=data.frame( cdf_hr7,x,z1,Z21,Z22)
tada8=tada7[order(tada7$cdf_hr7), ]
library('goftest')
F_sto=0
tt7=tada8$x
Z17=tada8$z1;
Z18=tada8$Z21;
Z19=tada8$Z22;

n=m
 MUu27=0
  for( e in 1:m)
  {
     fx7=function(z)
               {		
			z*(1-exp(-(z/(exp(B0Hr+B1Hr*Z17[e]+B21Hr*Z18[e]+B22Hr*Z19[e])))^(-AAHr)))						  
                  }
               tx7=integrate(fx7,lower=0,upper=w)
               MUu27[e]=tx7$value
  }            
 tx27=0
 for(h in 1:m)
 {
   fk13=function(z27)
                 {
			z27*(1-exp(-(z27/(exp(B0Hr+B1Hr*Z17[h]+B21Hr*Z18[h]+B22Hr*Z19[h])))^(-AAHr)))						  
                  }
   tx27[h]=integrate(fk13,lower=0,upper=tt7[h])$value
  }
  cdf_hr=tx27/MUu27
model_Hr = 1/(12*m) + sum((cdf_hr -((1:m)-.5)/m)^2)
CVM_Pval_Hr=1-pCvM(model_Hr, n = m)

##########################VREHT CvM
##############################################################################
######VREHT CvM######
 MUu2=0
  for( e in 1:m)
  {
     fx3=function(z)
               {					
			z*(1-(tanh(z/(exp(B0+B1*z1[e]+BB21*Z21[e]+BB22*Z22[e]))))^AA)						  
                  }
               tx81=integrate(fx3,lower=0,upper=w)
               MUu2[e]=tx81$value
  }            

 tx21=0
 for(h in 1:m)
 {
   fk11=function(z21)
                 {
			z21*(1-(tanh(z21/(exp(B0+B1*z1[h]+BB21*Z21[h]+BB22*Z22[h]))))^AA)						  
                  }
   tx21[h]=integrate(fk11,lower=0,upper=x[h])$value
  }
cdf_Tanho=tx21/MUu2
tada2=data.frame(cdf_Tanho,x,z1,Z21,Z22)
tada3=tada2[order(tada2$cdf_Tanho), ]
library('goftest')
F_sto=0
tt1=tada3$x
Zz13=tada3$z1;
Zz14=tada3$Z21;
Zz15=tada3$Z22;


######VREHT CvM######
 MUu4=0
  for( e in 1:m)
  {
     fx3=function(z)
               {					
			z*(1-(tanh(z/(exp(B0+B1*Zz13[e]+BB21*Zz14[e]+BB22*Zz15[e]))))^AA)						  
                  }
               tx=integrate(fx3,lower=0,upper=w)
               MUu4[e]=tx$value
  }            

 tx24=0
 for(h in 1:m)
 {
   fk14=function(z24)
                 {
			z24*(1-(tanh(z24/(exp(B0+B1*Zz13[h]+BB21*Zz14[h]+BB22*Zz15[h]))))^AA)						  
                  }
   tx24[h]=integrate(fk14,lower=0,upper=tt1[h])$value
  }
  cdf_Tanho4=tx24/MUu4
model_VREHT = 1/(12*m) + sum((cdf_Tanho4 -((1:m)-.5)/m)^2)
CVM_Pval_vreht=1-pCvM(model_VREHT, n = m)

############################################## Result##########

print(summary(resHn))
print(summary(resHr))
print(summary(res))

LLLLTn=res$maximum
LLLLhn=resHn$maximum
LLLHr=resHr$maximum
      
	  AIC_Hn=(-2*LLLLhn)+2*length(resHn$estimate);AIC_Hr=(-2*LLLHr)+2*length(resHr$estimate);AIC_TN=(-2*LLLLTn)+2*length(res$estimate);
	
	AIC_TN1[J]=(-2*LLLLTn)+2*length(res$estimate)
	AIC_Hn1[J]=(-2*LLLLhn)+2*length(resHn$estimate)
	AIC_Hr1[J]=(-2*LLLHr)+4*length(resHr$estimate)

	MM[J]=m
	D_Hn1[J]=D_Hn
	D_HR1[J]=D_HR
	D_vreht1[J]=D_vreht

      D_Hn19[J]=D_Hn18
	D_HR19[J]=D_HR18
	D_vreht19[J]=D_vreht18

      piHn3[J]=piHn12
	piHr3[J]=piHr12
	piVr3[J]=piTnL12

      pvalHn[J]=CVM_Pval_Hn;pvalHr[J]=CVM_Pval_Hr;pvalVRHT[J]=CVM_Pval_vreht;
      pvHN[J]=ifelse(CVM_Pval_Hn<=0.05,0,1)
      pvHR[J]=ifelse(CVM_Pval_Hr<=0.05,0,1)
      pvVHRT[J]=ifelse(CVM_Pval_vreht<=0.05,0,1)

    ########### to find delta AIC#####
     an1=AIC_Hn-min(AIC_Hn,AIC_Hr,AIC_TN)
     an2=AIC_Hr-min(AIC_Hn,AIC_Hr,AIC_TN)
     an3=AIC_TN-min(AIC_Hn,AIC_Hr,AIC_TN)

    AN6[J]=ifelse((an1<=2) && (CVM_Pval_Hn > 0.05) ,1,0);AN7[J]=ifelse((an2<=2) &&(CVM_Pval_Hr > 0.05),1,0);AN8[J]=ifelse(an3<=2 && (CVM_Pval_vreht> 0.05),1,0)     
      if((D_HR>=2.6) || (D_vreht>=2.6))
         {
           J=J
          }else{        
	J=J+1
       }
    }
	next
}
D=data.frame(MM,True_Pa,AIC_Hn1,AIC_Hr1,AIC_TN1,pvalHn,pvalHr,pvalVRHT,D_Hn1,D_HR1,D_vreht1,D_Hn19,D_HR19,D_vreht19,piHn3,piHr3,piVr3)
NN2=length(D$MM);d5=colMeans(x=D)
mean(D_Hn1)

D5=D[which(pvalHn>=0.05),c(1,3,6,9,12,15)]
D6=D[which(pvalHr>=0.05),c(1,4,7,10,13,16)]
D7=D[which(pvalVRHT>=0.05),c(1,5,8,11,14,17)]

nn1=length(D5$D_Hn1);nn2=length(D6$D_HR1);nn3=length(D7$D_vreht1)
DHN_av=mean(D5$D_Hn1);DHN_var=var(D5$D_Hn1);
DHR_av=mean(D6$D_HR1);DHR_var=var(D6$D_HR1);
DVR_av=mean(D7$D_vreht1);DVR_var=var(D7$D_vreht1)
DHN_av=mean(D5$piHn3);DHN_var=var(D5$piHn3);
DHR_av=mean(D6$piHr3);DHR_var=var(D6$piHr3);
DVR_av=mean(D7$piVr3);DVR_var=var(D7$piVr3)



CILHN=((mean(D5$D_Hn19))+(sqrt((var(D5$D_Hn19)))/(sqrt(nn1)))*abs(qnorm(0.05/2)))-((mean(D5$D_Hn19))-(sqrt((var(D5$D_Hn19)))/(sqrt(nn1)))*abs(qnorm(0.05/2)))
CILHR=((mean(D6$D_HR19))+(sqrt((var(D6$D_HR19)))/(sqrt(nn2)))*abs(qnorm(0.05/2)))-((mean(D6$D_HR19))-(sqrt((var(D6$D_HR19)))/(sqrt(nn2)))*abs(qnorm(0.05/2)))
CILVR=((mean(D7$D_vreht19))+(sqrt((var(D7$D_vreht19)))/(sqrt(nn3)))*abs(qnorm(0.05/2)))-((mean(D7$D_vreht19))-(sqrt((var(D7$D_vreht19)))/(sqrt(nn3)))*abs(qnorm(0.05/2)))

CILHNDa=(DHN_av+(sqrt(DHN_var)/(sqrt(nn1)))*abs(qnorm(0.05/2)))-(DHN_av-(sqrt(DHN_var)/(sqrt(nn1)))*abs(qnorm(0.05/2)))
CILHRDa=(DHR_av+(sqrt(DHR_var)/(sqrt(nn2)))*abs(qnorm(0.05/2)))-(DHR_av-(sqrt(DHR_var)/(sqrt(nn2)))*abs(qnorm(0.05/2)))
CILVRDa=(DVR_av+(sqrt(DVR_var)/(sqrt(nn3)))*abs(qnorm(0.05/2)))-(DVR_av-(sqrt(DVR_var)/(sqrt(nn3)))*abs(qnorm(0.05/2)))

CILHNpi=(mean(D5$piHn3)+(sqrt(var(D5$piHn3))/(sqrt(nn1)))*abs(qnorm(0.05/2)))-(mean(D5$piHn3)-(sqrt(var(D5$piHn3))/(sqrt(nn1)))*abs(qnorm(0.05/2)))
CILHRpi=(mean(D6$piHr3)+(sqrt(var(D6$piHr3))/(sqrt(nn2)))*abs(qnorm(0.05/2)))-(mean(D6$piHr3)-(sqrt(var(D6$piHr3))/(sqrt(nn2)))*abs(qnorm(0.05/2)))
CILVRpi=(mean(D7$piVr3)+(sqrt(var(D7$piVr3))/(sqrt(nn3)))*abs(qnorm(0.05/2)))-(mean(D7$piVr3)-(sqrt(var(D7$piVr3))/(sqrt(nn3)))*abs(qnorm(0.05/2)))

re_HN=c("sample size"=round(mean(D5$MM)),"Hn PMS"=sum(AN6)*100/NN2,"Hn PMF"=sum(pvHN)*100/NN2,"CvM HN"=mean(D5$pvalHn),"D_HN"=mean(D5$D_Hn1),"D_HN_Bias"=mean(D5$D_Hn1-1),"MSE_D_HN"=((mean(D5$D_Hn1)-1)^2)+var(D5$D_Hn1),"CILHN"=CILHN,"D_HNDa"=mean(D5$D_Hn19),"D_HN_BiasDa"=mean(D5$D_Hn19-1),"MSE_D_HNDa"=((mean(D5$D_Hn19)-1)^2)+var(D5$D_Hn19),"CILHNDa"=CILHNDa,"Pi_HN"=mean(D5$piHn3),"pi_HN_Bias"=mean(D5$piHn3-True_Pa),"MSE_pi_HN"=((mean(D5$piHn3)-True_Pa)^2)+var(D5$piHn3),"CILHNpi"=CILHNpi)
re_HR=c("sample size"=round(mean(D6$MM)),"Hr PMS"=sum(AN7)*100/NN2,"Hr PMF"=sum(pvHR)*100/NN2,"CvM HR"=mean(D6$pvalHr),"D_HR"=mean(D6$D_HR1),"D_HR_Bias"=mean(D6$D_HR1-1),"MSE_D_HR"=((mean(D6$D_HR1)-1)^2)+var(D6$D_HR1),"CILHR"=CILHR,"D_HRDa"=mean(D6$D_HR19),"D_HR_BiasDa"=mean(D6$D_HR19-1),"MSE_D_HRDa"=((mean(D6$D_HR19)-1)^2)+var(D6$D_HR19),"CILHRDa"=CILHRDa,"Pi_HR"=mean(D6$piHr3),"pi_HR_Bias"=mean(D6$piHr3-True_Pa),"MSE_pi_HR"=((mean(D6$piHr3)-True_Pa)^2)+var(D6$piHr3),"CILHRpi"=CILHRpi)
re_VR=c("sample size"=round(mean(D7$MM)),"VREHT PMS"=sum(AN8)*100/NN2,"VREHT PMF"=sum(pvVHRT)*100/NN2,"CvM VR"=mean(D7$pvalVRHT),"D_VR"=mean(D7$D_vreht1),"D_VR_Bias"=mean(D7$D_vreht1-1),"MSE_D_VR"=((mean(D7$D_vreht1)-1)^2)+var(D7$D_vreht1),"CILVR"=CILVR,"D_VRDa"=mean(D7$D_vreht19),"D_VR_BiasDa"=mean(D7$D_vreht19-1),"MSE_D_VRDa"=((mean(D7$D_vreht19)-1)^2)+var(D7$D_vreht19),"CILVRDa"=CILVRDa,"Pi_VR"=mean(D7$piVr3),"pi_VR_Bias"=mean(D7$piVr3-True_Pa),"MSE_pi_VR"=((mean(D7$piVr3)-True_Pa)^2)+var(D7$piVr3),"CILVRpi"=CILVRpi)

name=c("sample size","PMS","PMF","CvM_avg","D mean","D Bias","D MSE","D CIL","Da mean","Da Bias","Da MSE","Da CIL","Pa mean","Pa Bias","Pa MSE","Pa CIL")
g1=data.frame(name,re_HN,re_HR,re_VR,True_Pa)
g2=round(g1[,-1],digits=3)
NN2;

write.csv(D,  "results/MCDS_full_results.csv", row.names = FALSE)
write.csv(g2, "results/MCDS_summary_results.csv", row.names = FALSE)

