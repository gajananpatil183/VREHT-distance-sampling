# ============================================================
# CDS Simulation (Point-transect sampling): HN, HR and VREHT Detection Functions
# Author: Gajanan Patil and Shashibhushan Mhadik
# Description:
# This script performs a Monte Carlo simulation study under the
# Conventional Distance Sampling (CDS) framework for point transect data.
# Radial distance data are generated using a mixture detection function.
# The performance of Half-Normal (HN), Hazard-Rate (HR) and VREHT models
# is evaluated using AIC, Cramér–von Mises goodness-of-fit test,
# detection probability (Pa), confidence intervals and coverage probability.
# ============================================================

rm(list=ls())
# ---------------------------
# Required Libraries
# ---------------------------
library(mrds)
library(maxLik)
library(Distance)
library(Deriv)
library(goftest)
library(nleqslv)

# ---------------------------
# Parameters
# ---------------------------
aa = 25   # HR scale parameter
bb = 3    # HR shape parameter
si = 5    # HN scale parameter

n = 40    # Sample size
w = 50    # Truncation distance
nn=501    #Number of simulation

# Mixture model coefficients
Phi1 = 0.5 
Phi2 = 0.5

TNN1=0;TNN2=0;TNN3=0;code11=0
AN1=0;AN2=0;AN3=0;SLA=0;pvHN=0;pvVHRT=0;pvHR=0;True_Pa=0;HN_pa=0;HR_pa=0;VHRT_pa=0;HN_AIC=0;HR_AIC=0;VHRT_AIC=0;HN_Cvm=0;HR_Cvm=0;VHRT_Cvm=0
AN5=0;AN6=0;AN7=0;AN8=0;L_HN=0;U_HN=0;CL_HN=0;CP_HN=0;L_HR=0;U_HR=0;CL_HR=0;CP_HR=0;L_VHRT=0;U_VHRT=0;CL_VHRT=0;CP_VHRT=0

 E=1
  while(E<nn)
   {
     FF1=function(ZZ)  ######### to get the value of mu
         {
          ((Phi1)*(ZZ)*(exp(-(ZZ^2)/(2*si^2))))+((Phi2)*(ZZ)*(1-exp(-(ZZ/aa)^(-bb))))
          }
        	Vby2pii=integrate(FF1,lower=0,upper=w)$value

	distance=0 ###### to store values of distance sample data
      v2=1
	while(v2<n+1) 
 		{
   		repeat{
   		     u=runif(1,0,1) 
  		     gp=function(x){
       	        fn1=function(q)
         	         {
                       ((Phi1)*(q)*(exp(-(q^2)/(2*si^2))))+((Phi2)*(q)*(1-exp(-(q/aa)^(-bb))))
          	         }
          	Tt=tryCatch(integrate(fn1,lower=0,upper=x),error=function(Tt){})
		    if(length(Tt)==0)
 		         {
                      dd=w+5
 		         }else{
   	 	          Tt=Tt
  	               }
       	    dd1=Tt$value-u*Vby2pii
        	    return(dd1)
         	      }
            dd=tryCatch(nleqslv(5,gp)$x,error=function(dd){})
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
       x=distance
  #print(x)
  VVhat=2*pi*Vby2pii
  H0_hat=(2*pi)/VVhat
  True_pa=2/(w*w*H0_hat)



##################To analyze this data##########
 ############# HN  ##########
          fhn=0;
		qqhn=function(par)
     		{
       		Si=par[1]
	 		f2=function(z)
               	{
                   z*exp(-(z^2)/(2*Si^2))	  
                  }
               	ttHn=integrate(f2,lower=0,upper=w)
               	vby2pihn=ttHn$value
           	for( j in 1:length(x))
			{
            	fhn[j]=(x[j]*exp(-(x[j]^2)/(2*Si^2)))
            	ifelse(fhn[j]<1e-300,1e-300,fhn[j])        
			}
	    	wwhn=prod(fhn) 
          	ww1hn=ifelse(wwhn<1e-300,1e-300,wwhn)
		f3hn=(log(ww1hn)-length(x)*log(vby2pihn))   
          	return(f3hn)
     		}
		Ahn <- matrix(c(Si=1), 1, 1,byrow=TRUE)
		Bhn <- matrix(c(0),1,1)
		constraintshn <- list(ineqA=Ahn, ineqB=Bhn)
		reshn <- maxBFGS(qqhn, start=c(20), constraints=constraintshn)
		print(summary(reshn))
		sa1=reshn$estimate
		LLLhn=reshn$maximum
		si6=sa1[1]
		f32=function(z)
               {
                  z*exp(-(z^2)/(2*si6^2))  
                  }
               tt6=integrate(f32,lower=0,upper=w)
               v_hat6=2*pi*tt6$value
		   h0_hat6=(2*pi)/v_hat6
		   Pa_hat6=2/(w*w*h0_hat6)
	        AIC_pt_HN=(-2*LLLhn)+2


############# HR  ##########
      fhr=0;
		qqhr=function(par)
     		{
       	SI=par[1]
       	b=par[2]
	 	f2=function(z)
               {
                  (1-exp(-(z/SI)^(-b)))*z	  
                  }
               tthr=integrate(f2,lower=0,upper=w)
               vby2piHR=tthr$value
           	for( j in 1:length(x))
		{
            fhr[j]=((1-exp(-(x[j]/SI)^(-b)))*x[j])
            ifelse(fhr[j]<1e-300,1e-300,fhr[j])        
		}
	    	wwhr=prod(fhr) 
          	ww1hr=ifelse(wwhr<1e-300,1e-300,wwhr)
		f3hr=(log(ww1hr)-length(x)*log(vby2piHR))  
          	return(f3hr)
     		}
		A9 <- matrix(c(SI=1,b=0,SI=0,b=1), 2, 2,byrow=TRUE)
		B9 <- matrix(c(0,0),2,1,byrow=TRUE)
		constraintshr<- list(ineqA=A9, ineqB=B9)
		res9 <- maxBFGS(qqhr, start=c(20,5), constraints=constraintshr)
		print(summary(res9))
		sa3=res9$estimate
		LLLhr=res9$maximum
		SIHr=sa3[1]
		BHr=sa3[2]
			f3h=function(z)
               	{
                  	(1-exp(-(z/SIHr)^(-BHr)))*z		  
                  	}	
               	tt3hr=integrate(f3h,lower=0,upper=w)
             v_hathr=2*pi*tt3hr$value
		h0_hathr=(2*pi)/v_hathr
		Pa_hathr=2/(w*w*h0_hathr)
		Row_hathr=w*sqrt(Pa_hathr)
		AIC_pt_HR=(-2*LLLhr)+4
	
############### VHRT###############
fvh=0;
	qq=function(par)
     		{
       	 s=par[1]
       	 a=par[2]
	 	  f2vh=function(z)
                {
                  (1-(tanh(z/s))^a)*z	  
                  }
                 ttvh=integrate(f2vh,lower=0,upper=w)
                 vby2pivh=ttvh$value
              for( j in 1:length(x))
		   {
                fvh[j]=(x[j]*(1-(tanh(x[j]/s))^a))
                ifelse(fvh[j]<1e-300,1e-300,fvh[j])        
		   }
	          ww=prod(fvh) 
                ww1=ifelse(ww<1e-300,1e-300,ww)
		    f3=(log(ww1)-length(x)*log(vby2pivh))   
            return(f3)
          }
	Avh <- matrix(c(s=1,a=0,s=0,a=1), 2, 2,byrow=TRUE)
	Bvh <- matrix(c(0,-2),2,1)
	constraintsvh <- list(ineqA=Avh, ineqB=Bvh)
	res <- maxBFGS(qq, start=c(20,5), constraints=constraintsvh)
	print(summary(res))
	savh=res$estimate
	LLLvh=res$maximum
	Svh=savh[1]
	Avh=savh[2]
	f3=function(z)
               {
                  (1-(tanh(z/Svh))^Avh)*z	  
                  }
               ttvh2=integrate(f3,lower=0,upper=w)
               v_hat=2*pi*ttvh2$value
	h0_hat=(2*pi)/v_hat
	Pa_hat=2/(w*w*h0_hat)
	Row_hat=w*sqrt(Pa_hat)
	AIC_pt_Tn=(-2*LLLvh)+4


 ############# CvM Test P_Values
  	######  CvM test HN plot in point transect sampling##########
	library('goftest')
	F_sto2=0
	xb92=x
	yy102=sort(xb92)
	m112=length(yy102)
	n=m112
	######Origional KS_D######
  	y2=sort(x)
   	fo22=function(z)
               {
                  z*(exp(-(z^2)/(2*si6^2)))	  
                  }
	tt1122=integrate(fo22,lower=0,upper=w)
	ho2=tt1122$value
	tt122=0
 	for(h in 1:n)
 	{
   	fk12=function(z2)
               {
                  z2*(exp(-(z2^2)/(2*si6^2)))	  
                  }
  	 ttt2=integrate(fk12,lower=0,upper=y2[h])
   	tt122[h]=ttt2$value
  		}
  		cdf_hn=tt122/ho2
 	ghg2=0
	for(ss in 1:n)
   	{
     	ghg2[ss]=(2*ss-1)/(2*n)
   	}
	F_empo2=ghg2
	 D_ksO2=max(abs(F_empo2-cdf_hn))
 	 D_CMVhn=sum((cdf_hn-F_empo2)^2)
  	CVM_Pvalhn=1-pCvM(D_CMVhn, n =m112)
  	

###### KS and CvM test, in line transect sampling HR##########
      F_sto=0
	xb93=x
	yy103=sort(xb93)
	m11=length(yy103)
	n=m11
	#####empirical cdf###
 	ghg2=0
		for(ss in 1:m11)
   		{
    		 ghg2[ss]=(2*ss-1)/(2*n)
   		}
	F_empo5=ghg2

  ######Origional CvM######
  		y333=sort(x)
            y3=ifelse(y333==0,0.00000001,y333)
             so3=SIHr;bo=BHr
		fo23=function(z)
               {
                z*(1-exp(-(z/so3)^(-bo)))	  
                }
		tt1123=integrate(fo23,lower=0,upper=w)
		ho3=tt1123$value
		tt123=0
 	for(h in 1:n)
 		{
   		fk13=function(z2)
             {
             z2*(1-exp(-(z2/so3)^(-bo)))  
             }
   		tt3t=integrate(fk13,lower=0,upper=y3[h])
   	tt123[h]=tt3t$value
  	}
  	cdf_hr=tt123/ho3
	D_ksOhr=max(abs(F_empo5-cdf_hr))
	CvMhr=(1/(12*n))+sum((cdf_hr-F_empo5)^2)
	CVM_Pvalhr=1-pCvM(CvMhr, n = n) 


###### KS and CvM test,PP plot in line transect sampling Tan##########
  F_st=0
	xb91=x
	yy91=sort(xb91)
	m12=length(yy91)
	n=m12
	#####empirical cdf###
 	 gk1=0
		for(ss in 1:m12)
   		{
    		 gk1[ss]=(2*ss-1)/(2*n)
   		}
	F_empi=gk1

  ######Origional CvM######
  		y33=sort(x)
            So=Svh;Ao=Avh
		fo2=function(z)
               {
                z*(1-(tanh(z/So))^Ao)	  
                }
		tt5=integrate(fo2,lower=0,upper=w)
		ho5=tt5$value
		tt12=0
 	for(h in 1:n)
 		{
   		fk2=function(z22)
             {
             z22*(1-(tanh(z22/So))^Ao)  
             }
   		tt7=integrate(fk2,lower=0,upper=y33[h])
   	tt12[h]=tt7$value
  	}
  	cdf_TN=tt12/ho5
      D_ksOTn=max(abs(F_empi-cdf_TN))
      CvMTn=(1/(12*n))+sum((cdf_TN-F_empi)^2)
      CVM_PvalTn=1-pCvM(CvMTn, n = n)
      


############# To calaculate CI for Pa#############
###################### CI for Pa HN############### 
            r=x;si66=si6^2;n2=length(r)
		f51=function(z)
               	{
                 (exp(-(z^2)/(2*si66)))*z	  
                  }
               	AABB=integrate(f51,lower=0,upper=w)$value
		#####w r t si#####
		f52=function(si66)
           		{
            	log((exp(-(r^2)/(2*si66)))*r)	  
           		}
   		v52=Deriv(f52)
		f53=function(si66)
           		{
            	(exp(-(r^2)/(2*si66)))*r  
           		}
   		v53=Deriv(f53)
		f54=function(r)
        		{
        		 2 * (r^3 * exp(-(r^2/(2 * si66)))/(2 * si66)^2)
        	 	}
    		v54=integrate(f54,lower=0,upper=w)$value
		Dsi=(2 * (r^2/(2 * si66)^2))-(v54/AABB)
		######Hss
		Hsi1si1=(1/n2)*sum((Dsi)*(Dsi))
	##### Hessian is
	Hen5=matrix(c(Hsi1si1),ncol=1,byrow=TRUE)
	#######h0_s and h0_a
	h0_si=-v54/AABB^2
	h0_dash3=c(h0_si)
	Var_pt_hn=(1/n2)*t(h0_dash3)%*%solve(Hen5)%*%t(t(h0_dash3))
	se_hn_h0=sqrt(Var_pt_hn)
	CV_hn_h0=as.vector(se_hn_h0/(h0_hat6))

	###########
      CV_Denpt_hn=CV_hn_h0
	d.f.hn=(CV_Denpt_hn^4)/ as.vector((((CV_hn_h0^4)/(n2-1))))
	Chn=exp(qt(0.975,d.f.hn)*sqrt(log(1+((CV_hn_h0^2)))))
	Ll_hn_N=Pa_hat6/Chn
	Ul_hn_N=Pa_hat6*Chn



################## CI for Pa HR###########
############CI cal of HR
                n2=length(r)
				f41=function(z)
               		{
                 		(1-exp(-(z/SIHr)^(-BHr)))*z	  
                  	}
               		AAB=integrate(f41,lower=0,upper=w)$value

		#####w r t SI#####
				f42=function(SIHr)
           			{
            		log((1-exp(-(r/SIHr)^(-BHr)))*r)	  
           			}
   				v42=Deriv(f42)
				f43=function(SIHr)
           			{
            		(1-exp(-(r/SIHr)^(-BHr)))*r	  
           			}
  				 v43=Deriv(f43)
				f44=function(r)
        			{        
           			((BHr * r^2 * exp(-(r/SIHr)^-BHr))/(SIHr^2 * (r/SIHr)^(1 +BHr)))
         			}
   				 v44=integrate(f44,lower=0,upper=w)$value
				DSI=BHr* r * (exp(-(r/SIHr)^-BHr))/(SIHr^2 * (1 - (exp(-(r/SIHr)^-BHr))) * (r/SIHr)^(1 +BHr))-(v44/AAB)
   
		####### w r t b
				f45=function(BHr)
           			{
           			 log((1-exp(-(r/SIHr)^(-BHr)))*r)	  
           			}
   				v45=Deriv(f45)
				f46=function(BHr)
           			{
            		(1-exp(-(r/SIHr)^(-BHr)))*r	  
           			}
   				v46=Deriv(f46)
				f47=function(r)
        			{         
           			(-(r * exp(-(r/SIHr)^(-BHr))) * (log(r) - log(SIHr))/(r/SIHr)^BHr)
         			}
    				v47=integrate(f47,lower=0,upper=w)$value

				DB=-((exp(-(r/SIHr)^-BHr))* (log(r) - log(SIHr))/((1 - (exp(-(r/SIHr)^-BHr))) * (r/SIHr)^BHr))-(v47/AAB)
		######Hss
				Hsisi=(1/n2)*sum((DSI)*(DSI))
				HSIb=(1/n2)*sum((DSI)*(DB))
				HbSI=(1/n2)*sum((DB)*(DSI))
				Hbb=(1/n2)*sum((DB)*(DB))
		##### Hessian is
				HenHr=matrix(c(Hsisi,HSIb,HbSI,Hbb),ncol=2,byrow=TRUE)
				#######h0_s and h0_a
				h0_SI=-v44/AAB^2
				h0_b=-v47/AAB^2
				h0_dash2=c(h0_SI,h0_b)
				Var_pt_hr=(1/n2)*t(h0_dash2)%*%solve(HenHr)%*%t(t(h0_dash2))
				se_hr_h0=sqrt(Var_pt_hr)
				CV_hr_h0=as.vector(se_hr_h0/(h0_hathr))

		###########
      CV_Denpt_hr=CV_hr_h0
	d.f.hr=(CV_Denpt_hr^4)/ as.vector((((CV_hr_h0^4)/(n2-1))))
	Chr=exp(qt(0.975,d.f.hr)*sqrt(log(1+((CV_hr_h0^2)))))
	Ll_hr_N=Pa_hathr/Chr
	Ul_hr_N=Pa_hathr*Chr

#################### CI for Pa VHRT ##############


 library(Deriv)
	s=Svh
	a=Avh
	x=replace(x,x==0,0.000000000001)
	r=x
	n2=length(r)   
	f31=function(z)
               	{
                  (1-tanh(z/s)^a)*z	  
                  }
               	AAA=integrate(f31,lower=0,upper=w)$value
	#####w r t s#####
	f32=function(s)
           		{
            	log((1-tanh(r/s)^a)*r)	  
           		}
   			v32=Deriv(f32)
	f33=function(s)
           		{
           		((1-tanh(r/s)^a)*r)	  
           		}
   			v33=Deriv(f33)
	f34=function(r)
        		{
         		a * r^2 * (1 - (tanh(r/s))^2) *(tanh(r/s))^(a - 1)/s^2
         		}
   	 		v34=integrate(f34,lower=0,upper=w)$value
	DS=(a * r * (1 - (tanh(r/s))^2) * (tanh(r/s))^(a - 1)/(s^2 * (1 - (tanh(r/s))^a)))-(v34/AAA)
	########w r t a
	f35=function(a)
           		{
            	log((1-tanh(r/s)^a)*r)	  
           		}
   			v35=Deriv(f35)
	f36=function(a)
           		{
            	((1-tanh(r/s)^a)*r)	  
           		}
  			v36=Deriv(f36)
	f37=function(r)
        		{
        		-(r * log(tanh(r/s)) * (tanh(r/s))^a)
         		}
    			v37=integrate(f37,lower=0,upper=w)$value
	DA=-(log(tanh(r/s))*((tanh(r/s))^a)/(1 - ((tanh(r/s))^a)))-(v37/AAA)
	######Hss
		Hss=(1/n2)*sum((DS)*(DS))
		Hsa=(1/n2)*sum((DS)*(DA))
		Has=(1/n2)*sum((DA)*(DS))
		Haa=(1/n2)*sum((DA)*(DA))
	##### Hessian is
		Hessian1=matrix(c(Hss,Hsa,Has,Haa),ncol=2,byrow=TRUE)
	#######h0_s and h0_a
		h0_s=-v34/AAA^2
		h0_a=-v37/AAA^2
		h0_dash=c(h0_s,h0_a)
		Var_pt_tn=(1/n2)*t(h0_dash)%*%solve(Hessian1)%*%t(t(h0_dash))
		se_h0=sqrt(Var_pt_tn)
		CV_vrht_h0=se_h0/(h0_hat)

	###########
      CV_Denpt_vrht=CV_vrht_h0
	d.f.vr=(CV_Denpt_vrht^4)/ as.vector((((CV_vrht_h0^4)/(n2-1))))
	Chvr=exp(qt(0.975,d.f.vr)*sqrt(log(1+((CV_vrht_h0^2)))))
	Ll_vrht_N=Pa_hat/Chvr
	Ul_vrht_N=Pa_hat*Chvr


        Code1=res$code;Code2=res9$code
        
	   if(Code1==0  && Code2==0 )
        {

        code11[E]=Code1
#############################
      True_Pa=True_pa
      HN_pa[E]=Pa_hat6
      HR_pa[E]=Pa_hathr
      VHRT_pa[E]=Pa_hat

      HN_AIC[E]=AIC_pt_HN
      HR_AIC[E]=AIC_pt_HR
      VHRT_AIC[E]=AIC_pt_Tn
  
      HN_Cvm[E]=CVM_Pvalhn
      HR_Cvm[E]=CVM_Pvalhr
      VHRT_Cvm[E]=CVM_PvalTn


      TNN1[E]=CV_hn_h0
      TNN2[E]=CV_hr_h0
      TNN3[E]=CV_vrht_h0
     
      pvHN[E]=ifelse(CVM_Pvalhn<=0.05,0,1)
      pvHR[E]=ifelse(CVM_Pvalhr<=0.05,0,1)
      pvVHRT[E]=ifelse(CVM_PvalTn<=0.05,0,1)

      L_HN[E]=Ll_hn_N
      U_HN[E]=Ul_hn_N
      CL_HN[E]=ifelse(Ul_hn_N>1,1,Ul_hn_N)-ifelse(Ll_hn_N<0,0,Ll_hn_N)
      CP_HN[E]=ifelse(Ll_hn_N <= True_Pa && True_Pa <= Ul_hn_N,1,0)

      L_HR[E]=Ll_hr_N
      U_HR[E]=Ul_hr_N
      CL_HR[E]=ifelse(Ul_hr_N>1,1,Ul_hr_N)-ifelse(Ll_hr_N<0,0,Ll_hr_N)
      CP_HR[E]=ifelse(Ll_hr_N <= True_Pa&& True_Pa <=Ul_hr_N,1,0)

      L_VHRT[E]=Ll_vrht_N
      U_VHRT[E]=Ul_vrht_N
      CL_VHRT[E]=ifelse(Ul_vrht_N>1,1,Ul_vrht_N)-ifelse(Ll_vrht_N<0,0,Ll_vrht_N)
      CP_VHRT[E]=ifelse(Ll_vrht_N <= True_Pa&& True_Pa <=Ul_vrht_N,1,0)
  

      #AN1[E]=ifelse(((-0.01 >= Res$gradient[1] || Res$gradient[1] >= 0.01) || (-0.01 >= Res$gradient[2] || Res$gradient[2] >= 0.01 ) ),1,0)
      #AN2[E]=ifelse((-0.01 >= Res$gradient[1] || Res$gradient[1] >= 0.01),1,0)
      #AN3[E]=ifelse((-0.01 >= Res$gradient[2] || Res$gradient[2] >= 0.01),1,0)
      #AN5[E]=Res$gradient[1];#AN6[E]=Res$gradient[2]

 
    ########### to find delta AIC#####
     an1=AIC_pt_HN-min(AIC_pt_HN,AIC_pt_HR,AIC_pt_Tn)
     an2=AIC_pt_HR-min(AIC_pt_HN,AIC_pt_HR,AIC_pt_Tn)
     an3=AIC_pt_Tn-min(AIC_pt_HN,AIC_pt_HR,AIC_pt_Tn)
     
    AN6[E]=ifelse((an1<=2) && (CVM_Pvalhn > 0.05) ,1,0);AN7[E]=ifelse((an2<=2) &&(CVM_Pvalhr > 0.05),1,0);AN8[E]=ifelse(an3<=2 && (CVM_PvalTn> 0.05),1,0)     

      p1hn=CVM_Pvalhn;p1hr=CVM_Pvalhr;p1tn=CVM_PvalTn;  
      AIC_Hn=AIC_pt_HN;AIC_Hr=AIC_pt_HR;AIC_TN =AIC_pt_Tn

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
      
	  next
      }	     	   
   }

mg=length(SLA)
ut=data.frame(SLA,True_Pa,HN_pa,HR_pa,VHRT_pa,HN_AIC,HR_AIC,VHRT_AIC,HN_Cvm,HR_Cvm,VHRT_Cvm,L_HN,U_HN,L_HR,U_HR,L_VHRT,U_VHRT,CL_HN,CL_HR,CL_VHRT,CP_HN,CP_HR,CP_VHRT,n,mg,TNN1,TNN2,TNN3)
sum(AN1);sum(AN2);sum(AN3);AN4=data.frame(AN5,AN6)
ut2=ut[order(SLA),]
gg5=colMeans(ut2)
MShn=length(which(SLA==1))
MShr=length(which(SLA==2))
MSvhrt=length(which(SLA==3))


ut88=data.frame(n,mg,"TRUE_Pa"=mean(True_Pa),"pahn"=mean(HN_pa),"pahr"=mean(HR_pa),"paVHRT"=mean(VHRT_pa),"MES_HN"=(((mean(HN_pa)-mean(True_Pa))^2)+var(HN_pa)),"MSE_HR"=(((mean(HR_pa)-mean(True_Pa))^2)+var(HR_pa)),"MSE_VHRT"=(((mean(VHRT_pa)-mean(True_Pa))^2)+var(VHRT_pa)),"HN fit"=(sum(pvHN)*100/mg),"HR fit"=(sum(pvHR)*100/mg),"VHRT fit"=(sum(pvVHRT)*100/mg),"CL HN"=mean(CL_HN),"CL HR"=mean(CL_HR),"CL VHRT"=mean(CL_VHRT),"CP HN"=sum(CP_HN)*100/mg,"CP HR"=sum(CP_HR)*100/mg,"CP VHRT"=sum(CP_VHRT)*100/mg)
ut88
mg
Report_HN=data.frame(True_Pa,HN_pa,L_HN,U_HN,CL_HN,CP_HN)
Report_HR=data.frame(True_Pa,HR_pa,L_HR,U_HR,CL_HR,CP_HR)
Report_VRHT=data.frame(True_Pa,VHRT_pa,L_VHRT,U_VHRT,CL_VHRT,CP_VHRT)
DaicHN=mean(AN6)*100;DaicHR=mean(AN7)*100;DaicVRHT=mean(AN8)*100;

ut71=ut[which(HN_Cvm>=0.05),c(1,9,3,12,13,18,21,26)]
biashn=(ut71$HN_pa-ut88$TRUE_Pa);mse_hn1=((ut71$TNN1*ut71$HN_pa)^2)+(biashn)^2
HN_rep=c(True_Pa,length(ut71[,1]),((length(ut71[,1])/mg)*100),colMeans(ut71[,c(-1,-7,-8)]),mean(biashn),mean(mse_hn1),"MES_Hn"=(((mean(ut71[,3])-mean(True_Pa))^2)+var(ut71[,3])),"CP_HN"=mean(ut71[,7])*100,DaicHN)

ut72=ut[which(HR_Cvm>=0.05),c(1,10,4,14,15,19,22,27)]
biashr=(ut72$HR_pa-ut88$TRUE_Pa);mse_hr1=((ut72$TNN2*ut72$HR_pa)^2)+(biashr)^2
HR_rep=c(True_Pa,length(ut72[,1]),((length(ut72[,1])/mg)*100),colMeans(ut72[,c(-1,-7,-8)]),mean(biashr),mean(mse_hr1),"MES_Hr"=(((mean(ut72[,3])-mean(True_Pa))^2)+var(ut72[,3])),"CP_HR"=mean(ut72[,7])*100,DaicHR)

ut73=ut[which(VHRT_Cvm>=0.05),c(1,11,5,16,17,20,23,28)]
biasvrht=(ut73$VHRT_pa-ut88$TRUE_Pa);mse_vrht1=((ut73$TNN3*ut73$VHRT_pa)^2)+(biasvrht)^2
VHRT_rep=c(True_Pa,length(ut73[,1]),((length(ut73[,1])/mg)*100),colMeans(ut73[,c(-1,-7,-8)]),mean(biasvrht),mean(mse_vrht1),"MES_Vhrt"=(((mean(ut73[,3])-mean(True_Pa))^2)+var(ut73[,3])),"CP_VHRT"=mean(ut73[,7])*100,DaicVRHT)

ut74=c("True Pa","total sample","perc. fit","p-value CvM","Pa_hat","ave LL","ave UL","ave CL","ave Bias","ave New MSE","ave MSE","CP","D AIC")

ut88
ut75=data.frame(ut74,HN_rep,HR_rep,VHRT_rep)
ut75

write.csv(ut, "results/full_simulation_results.csv", row.names = FALSE)
write.csv(ut88, "results/summary_results.csv", row.names = FALSE)
write.csv(ut75, "results/model_comparison.csv", row.names = FALSE)

