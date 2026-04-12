# ============================================================
# CDS Simulation(Line-transect sampling): HN, HR and VREHT Detection Functions
# Author: Gajanan Patil and Shashibhushan Mhadik
# Description:
# This script performs a Monte Carlo simulation study under the
# Conventional Distance Sampling (CDS) framework for line transect data
# comparing Half-Normal (HN), Hazard-Rate (HR) and VREHT detection functions
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

# ---------------------------
# Storage Initialization
# ---------------------------
TNN1=0;TNN2=0;TNN3=0
AN1=0;AN2=0;AN3=0;SLA=0;pvHN=0;pvVHRT=0;pvHR=0;True_Pa=0;HN_pa=0;HR_pa=0;VHRT_pa=0;HN_AIC=0;HR_AIC=0;VHRT_AIC=0;HN_Cvm=0;HR_Cvm=0;VHRT_Cvm=0
AN5=0;AN6=0;AN7=0;AN8=0;L_HN=0;U_HN=0;CL_HN=0;CP_HN=0;L_HR=0;U_HR=0;CL_HR=0;CP_HR=0;L_VHRT=0;U_VHRT=0;CL_VHRT=0;CP_VHRT=0

 E=1
  while(E<nn)
   {
     FF1=function(ZZ)  ######### to get the value of mu
         {
          ((Phi1)*(exp(-(ZZ^2)/(2*si^2))))+((Phi2)*(1-exp(-(ZZ/aa)^(-bb))))
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
                       ((Phi1)*(exp(-(q^2)/(2*si^2))))+((Phi2)*(1-exp(-(q/aa)^(-bb))))
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
       x=distance
  #print(x)
##################To analyze this data##########
 ############# HN  ##########
      f=0;
	qq2=function(par)
     {
       Si=par[1]
	 f22=function(z)
               	{
                 	exp(-(z^2)/(2*Si^2))	  
                  }
               	tt2=integrate(f22,lower=0,upper=w)
               	mu2=tt2$value
           for( j in 1:length(x))
			{
            	f[j]=(exp(-(x[j]^2)/(2*Si^2)))
            	ifelse(f[j]<1e-300,1e-300,f[j])        
			}
	    	ww2=prod(f) 
          	ww1=ifelse(ww2<1e-300,1e-300,ww2)
			f32=(log(ww1)-length(x)*log(mu2)) 
          		return(f32)
     	}
	Aa2 <- matrix(c(Si=1), 1, 1,byrow=TRUE)
	Bb2 <- matrix(c(0),1,1)
	constraintsHN <- list(ineqA=Aa2, ineqB=Bb2)
	res <- maxBFGS(qq2, start=c(20), constraints=constraintsHN)
	print(summary(res))
	sa=res$estimate
	LLL=res$maximum
	siHN=sa[1]
	f32=function(z)
               {
                  exp(-(z^2)/(2*siHN^2)) 	  
                  }
               tt2=integrate(f32,lower=0,upper=w)
               muN2=tt2$value
	f0_Hn=(1/muN2)
	mu_Hn=(1/f0_Hn)
	pa_Hn=mu_Hn/w
	AIC_Hn=(-2*LLL)+2


############# HR  ##########
fHR=0;
qq=function(par)
     {
       SI=par[1]
       bH=par[2]
	 f2=function(z)
               	{
                  (1-exp(-(z/SI)^(-bH)))	  
                  }
               	tthr1=tryCatch(integrate(f2,lower=0,upper=w),error=function(tthr1){})
               	mu17=tthr1$value
                  muTn18=tryCatch(log(mu17),error=function(mu17){})
                    if(length(muTn18)==0)
                        {
                            muN3=0                                                      
                      }else{
                         muN3=muTn18
                        }
           	for( j in 1:length(x))
			{
            	fHR[j]=((1-exp(-(x[j]/SI)^(-bH))))
            	ifelse(fHR[j]<1e-300,1e-300,fHR[j])        
			}
	    		wwhr1=prod(fHR) 
          		ww1hr1=ifelse(wwhr1<1e-300,1e-300,wwhr1)
			f3hr1=(log(ww1hr1)-length(x)*(muN3))   
          		return(f3hr1)
     			}
	A9 <- matrix(c(SI=1,b=0,SI=0,b=1), 2, 2,byrow=TRUE)
	B9 <- matrix(c(0,0),2,1,byrow=TRUE)
	constraintsHR <- list(ineqA=A9, ineqB=B9)
	res92<- maxBFGS(qq, start=c(20,5), constraints=constraintsHR)
      res91=tryCatch(res92,error=function(res92){})

	
############### VHRT###############
f2=0;
QQ=function(par)
     {
       s=par[1]
       a=par[2]
	 ff=function(z)
               	{
                  (1-(tanh(z/s))^a)	  
                  }
               	ttnn=tryCatch(integrate(ff,lower=0,upper=w),error=function(ttnn){})
               	muTn15=ttnn$value
                  muTn16=tryCatch(log(muTn15),error=function(muTn15){})
                    if(length(muTn16)==0)
                        {
                            muTn=0
                      }else{
                         muTn=muTn16
                        }
           	for( j in 1:length(x))
			{
            	f2[j]=(1-(tanh(x[j]/s))^a)
            	ifelse(f2[j]<1e-300,1e-300,f2[j])        
			}
	    		wwTn1=prod(f2) 
          		#wwTn=ifelse(wwTn1<1e-300,1e-300,wwTn1)
			V1=(log(wwTn1)-(length(x)*muTn))   
          		return(V1)
     			}
	Aa <- matrix(c(s=1,a=0,s=0,a=1), 2, 2,byrow=TRUE)
	Bb <- matrix(c(0,0),2,1,byrow=TRUE)
	constraintsab <- list(ineqA=Aa, ineqB=Bb)
	ess <- maxBFGS(QQ, start=c(20,3), constraints=constraintsab)
      Res=tryCatch(ess,error=function(ess){})
      if((length(Res)==0) || (length(res91)==0)|| (muN3=0) || (muTn=0)|| is.na(res91$gradient[1]) || is.na(res91$gradient[2]) || is.na(Res$gradient[1]) || is.na(Res$gradient[2]) || (-0.01 >= Res$gradient[1] || Res$gradient[1] >= 0.01) || (-0.01 >= Res$gradient[2] || Res$gradient[2] >= 0.01 ) || (-0.01 >= res91$gradient[1] || res91$gradient[1] >= 0.01) || (-0.01 >= res91$gradient[2] || res91$gradient[2] >= 0.01 ) )
      {
          E=E
       }else{
	print(summary(Res))
      print(summary(res91))
	sa3=res91$estimate
	LLLr=res91$maximum
	SIHR=sa3[1]
	bHR=sa3[2]
       fm3=function(z)
               {
                  (1-exp(-(z/SIHR)^(-bHR)))		  
                  }
               tt33=integrate(fm3,lower=0,upper=w)
               mu_hathr=tt33$value
	f0_Hr=(1/mu_hathr)
	mu_Hr=(1/f0_Hr)
	pa_Hr=mu_Hr/w
	AIC_Hr=(-2*LLLr)+4

      #Res <- maxLik(QQ, start=c(20,5), constraints=constraintsab)
	#print(summary(Res))
	eP=Res$estimate
	Lval=Res$maximum
	Svh=eP[1]
	Avh=eP[2]
	ff2=function(z)
               {
                  1-(tanh(z/Svh))^Avh		  
                  }
               ttn2=integrate(ff2,lower=0,upper=w)
               muTN=ttn2$value
	f0_Tn=(1/muTN)
	mu_TN=(1/f0_Tn)
	pa_TN=mu_TN/w
	AIC_TN=(-2*Lval)+4


############# CvM Test P_Values
  ###### KS and CvM test,PP plot in line transect sampling##########
  library('goftest')
  F_sto2=0
  xbb=x
  yhln=sort(xbb)
  m112=length(yhln)
  n=m112
	######Origional KS_D######
  			y2=sort(x)
   			flhn=function(z)
               		{
                  	(exp(-(z^2)/(2*siHN^2)))	  
                  	}
			tlhn=integrate(flhn,lower=0,upper=w)
			hlhn=tlhn$value
			tlhn22=0
 		for(h in 1:n)
 			{
   			flhn2=function(z2)
               	{
                  (exp(-(z2^2)/(2*siHN^2)))	  
                  }
   			tlhn2=integrate(flhn2,lower=0,upper=y2[h])
   			tlhn22[h]=tlhn2$value
  			}
  			cdf_hnL=tlhn22/hlhn
 		ghln=0
		for(ss in 1:n)
   			{
     			ghln[ss]=(2*ss-1)/(2*n)
   			}
		F_empo2L=ghln
		D_ksO2L=max(abs(F_empo2L-cdf_hnL))
		D_CMVhnL=sum((cdf_hnL-F_empo2L)^2)
		CVM_PvalhnL=1-pCvM(D_CMVhnL, n = n)

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
             so3=SIHR;bo=bHR
		fo23=function(z)
               {
                (1-exp(-(z/so3)^(-bo)))	  
                }
		tt1123=integrate(fo23,lower=0,upper=w)
		ho3=tt1123$value
		tt123=0
 	for(h in 1:n)
 		{
   		fk13=function(z2)
             {
             (1-exp(-(z2/so3)^(-bo)))  
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
                (1-(tanh(z/So))^Ao)	  
                }
		tt5=integrate(fo2,lower=0,upper=w)
		ho5=tt5$value
		tt12=0
 	for(h in 1:n)
 		{
   		fk2=function(z2)
             {
             (1-(tanh(z2/So))^Ao)  
             }
   		tt7=integrate(fk2,lower=0,upper=y33[h])
   	tt12[h]=tt7$value
  	}
  	cdf_TN=tt12/ho5
      D_ksOTn=max(abs(F_empi-cdf_TN))
      CvMTn=(1/(12*n))+sum((cdf_TN-F_empi)^2)
      CVM_PvalTn=1-pCvM(CvMTn, n = n)
      CVM_PvalTn


############# To calaculate CI for Pa#############
###################### CI for Pa HN############### 
            r=x
		n2=length(x)
            n22=n2
		f51=function(z)
               	{
                 exp(-(z^2)/(2*siHN^2)) 	  
                  }
               	AABB=integrate(f51,lower=0,upper=w)$value
		#####w r t si#####
		f52=function(siHN)
           		{
            	log(exp(-(z^2)/(2*siHN^2)) )	  
           		}
   		v52=Deriv(f52)
		f53=function(siHN)
           		{
            	(exp(-(z^2)/(2*siHN^2))) 
           		}
   		v53=Deriv(f53)
		f54=function(z)
        		{ 
    			4 * (siHN * (z^2) * exp(-((z^2)/(2 * siHN^2)))/(2 * siHN^2)^2)
        	 	}
    		v54=integrate(f54,lower=0,upper=w)$value
		Dsi=(4 * (siHN * x^2/(2 * siHN^2)^2))-(v54/AABB)
		######Hss
		Hsi1si1=(1/n22)*sum((Dsi)*(Dsi))
	##### Hessian is
	Hen5=matrix(c(Hsi1si1),ncol=1,byrow=TRUE)
	#######h0_s and h0_a
	f0_si=-v54/AABB^2
	f0_dash3=c(f0_si)
	Var_lt_hn=(1/n2)*t(f0_dash3)%*%solve(Hen5)%*%t(t(f0_dash3))
	se_hn_f0=sqrt(Var_lt_hn)
	CV_hn_h0=as.vector(se_hn_f0/(f0_Hn))
	#####Var Hn#########
	CV_Den=CV_hn_h0
	d.f.=CV_Den^4 / ((CV_hn_h0^4/(n22-1)))
	C=exp(qt(0.975,d.f.)*sqrt(log(1+((CV_hn_h0^2)))))
	Llhn=pa_Hn/C
	Ulhn=pa_Hn*C


################## CI for Pa HR###########
############CI cal of HR

		 n255=length(r)
		 f41=function(z)
               {
                 (1-exp(-(z/SIHR)^(-bHR)))	  
                  }
               AAB=integrate(f41,lower=0,upper=w)$value

        #####w r t SI#####
        f42=function(SIHR)
           	{
            log((1-exp(-(r/SIHR)^(-bHR))))	  
           	}
   		v42=Deriv(f42)
		f43=function(SIHR)
          	{
            (1-exp(-(r/SIHR)^(-bHR)))	  
           	}
   		v43=Deriv(f43)
		f44=function(r)
        	{
           	( bHR * r * exp(-(r/SIHR)^-bHR)/(SIHR^2 * (r/SIHR)^(1 + bHR)))
         	}
         v44=integrate(f44,lower=0,upper=w)$value
	DSI=(bHR * r * (exp(-(r/SIHR)^-bHR))/(SIHR^2 * (1 - (exp(-(r/SIHR)^-bHR))) * (r/SIHR)^(1 + bHR)))-(v44/AAB)
	####### w r t b
	f45=function(bHR)
           	{
            log((1-exp(-(r/SIHR)^(-bHR))))	  
           	}
   		v45=Deriv(f45)
	f46=function(bHR)
           	{
            (1-exp(-(r/SIHR)^(-bHR)))	  
           	}
   		v46=Deriv(f46)

		f47=function(r)
        	{
           	( -(exp(-(r/SIHR)^-bHR) * (log(r) - log(SIHR))/(r/SIHR)^bHR) )
         	}
    		v47=integrate(f47,lower=0,upper=w)$value
		DB=(-((exp(-(r/SIHR)^-bHR)) * (log(r) - log(SIHR))/((1 - (exp(-(r/SIHR)^-bHR))) *(r/SIHR)^bHR)))-(v47/AAB)

		######Hss
			Hsisi=(1/n255)*sum((DSI)*(DSI))
			HSIb=(1/n255)*sum((DSI)*(DB))
			HbSI=(1/n255)*sum((DB)*(DSI))
			Hbb=(1/n255)*sum((DB)*(DB))
		##### Hessian is
		Hen=matrix(c(Hsisi,HSIb,HbSI,Hbb),ncol=2,byrow=TRUE)
		#######h0_s and h0_a
			f0_SI=-v44/AAB^2
			f0_b=-v47/AAB^2
			f0_dash2=c(f0_SI,f0_b)

		Var_lt_hr=(1/n255)*t(f0_dash2)%*%solve(Hen)%*%t(t(f0_dash2))
		se_hr_f0=sqrt(Var_lt_hr)
		CV_hr_f0=as.vector(se_hr_f0/(f0_Hr))

		#####Var Hn#########
			CV_Denhr=CV_hr_f0
			d.f.hr=CV_Denhr^4 /((CV_hr_f0^4/(n255-1)))
			Chr=exp(qt(0.975,d.f.hr)*sqrt(log(1+((CV_hr_f0^2)))))
		Llhr=pa_Hr/C
	      Ulhr=pa_Hr*C

#################### CI for Pa VHRT ##############


##########CI Tanh#######
		r=ifelse(x==0,0.0001,x)
		nn1=length(r)
		f5=function(z)
               	{
                 	(1-(tanh(z/Svh))^Avh)	  
                  }
               	AaBb=integrate(f5,lower=0,upper=w)$value

		#####w r t S#####
			f6=function(Svh)
           		{
            	log(1-(tanh(z/Svh))^Avh)	  
           		}
   			v6=Deriv(f6)
		f7=function(Svh)
          		{
            	((1-(tanh(z/Svh))^Avh))	  
           		}
   			v7=Deriv(f7)
		f8=function(r)
        		{
           		( Avh * r * (1 - (tanh(r/Svh))^2) *(tanh(r/Svh))^(Avh - 1)/Svh^2)
         		}
     		v8=integrate(f8,lower=0,upper=w)$value
		DS=(Avh * r * (1 - (tanh(r/Svh))^2) * (tanh(r/Svh))^(Avh - 1)/(Svh^2 * (1 - (tanh(r/Svh))^Avh)))-(v8/AaBb)
	####### w r t A########
		f9=function(Avh)
           	{
            log(1-(tanh(z/Svh))^Avh)	  
           	}
   		v9=Deriv(f9)
		f10=function(Avh)
           	{
            ((1-(tanh(z/Svh))^Avh))	  
           	}
   		v10=Deriv(f10)
		f11=function(r)
        	{
           	( -(log(tanh(r/Svh))*(tanh(r/Svh))^Avh))
         	}
    		v11=integrate(f11,lower=0,upper=w)$value
	DA=(-(log(tanh(r/Svh)) * ((tanh(r/Svh))^Avh)/(1 - ((tanh(r/Svh))^Avh))))-(v11/AaBb)

	######Hss
		HSS=(1/nn1)*sum((DS)*(DS))
		HSA=(1/nn1)*sum((DS)*(DA))
		HAS=(1/nn1)*sum((DA)*(DS))
		HAA=(1/nn1)*sum((DA)*(DA))
	##### Hessian is######
		HenTn=matrix(c(HSS,HSA,HAS,HAA),ncol=2,byrow=TRUE)
	#######f0_s and f0_a###
		f0_s=-v8/AaBb^2
		f0_a=-v11/AaBb^2
		f0_daSh=c(f0_s,f0_a)

	Var_lt_Tn=(1/nn1)*t(f0_daSh)%*%solve(HenTn)%*%t(t(f0_daSh))
	se_Tn_f0=sqrt(Var_lt_Tn)
	CV_Tn_f0=as.vector(se_Tn_f0/(f0_Tn))

	#####Var Tn#########
	CV_DenTN=CV_Tn_f0
	d.f.=CV_DenTN^4 / ((CV_Tn_f0^4/(nn1-1)))
	C=exp(qt(0.975,d.f.)*sqrt(log(1+((CV_Tn_f0^2)))))
	Llvhrt=pa_TN/C
	Ulvhrt=pa_TN*C



#############################
      True_Pa=MUU/w
      HN_pa[E]=pa_Hn
      HR_pa[E]=pa_Hr
      VHRT_pa[E]=pa_TN

      HN_AIC[E]=AIC_Hn
      HR_AIC[E]=AIC_Hr
      VHRT_AIC[E]=AIC_TN
  
      HN_Cvm[E]=CVM_PvalhnL
      HR_Cvm[E]=CVM_Pvalhr
      VHRT_Cvm[E]=CVM_PvalTn


      TNN1[E]=CV_hn_h0
      TNN2[E]=CV_hr_f0
      TNN3[E]=CV_Tn_f0
     
      pvHN[E]=ifelse(CVM_PvalhnL<=0.05,0,1)
      pvHR[E]=ifelse(CVM_Pvalhr<=0.05,0,1)
      pvVHRT[E]=ifelse(CVM_PvalTn<=0.05,0,1)

      L_HN[E]=Llhn
      U_HN[E]=Ulhn
      CL_HN[E]=Ulhn-Llhn
      CP_HN[E]=ifelse(Llhn <= True_Pa && True_Pa <= Ulhn,1,0)

      L_HR[E]=Llhr
      U_HR[E]=Ulhr
      CL_HR[E]=Ulhr-Llhr
      CP_HR[E]=ifelse(Llhr <= True_Pa&& True_Pa <=Ulhr,1,0)

      L_VHRT[E]=Llvhrt
      U_VHRT[E]=Ulvhrt
      CL_VHRT[E]=Ulvhrt-Llvhrt
      CP_VHRT[E]=ifelse(Llvhrt <= True_Pa&& True_Pa <=Ulvhrt,1,0)
  

      AN1[E]=ifelse(((-0.01 >= Res$gradient[1] || Res$gradient[1] >= 0.01) || (-0.01 >= Res$gradient[2] || Res$gradient[2] >= 0.01 ) ),1,0)
      AN2[E]=ifelse((-0.01 >= Res$gradient[1] || Res$gradient[1] >= 0.01),1,0)
      AN3[E]=ifelse((-0.01 >= Res$gradient[2] || Res$gradient[2] >= 0.01),1,0)
      AN5[E]=Res$gradient[1];AN6[E]=Res$gradient[2]

   
    ########### to find delta AIC#####
     an1=AIC_Hn-min(AIC_Hn,AIC_Hr,AIC_TN)
     an2=AIC_Hr-min(AIC_Hn,AIC_Hr,AIC_TN)
     an3=AIC_TN-min(AIC_Hn,AIC_Hr,AIC_TN)
     
    AN6[E]=ifelse((an1<=2) && (CVM_PvalhnL > 0.05) ,1,0);AN7[E]=ifelse((an2<=2) &&(CVM_Pvalhr > 0.05),1,0);AN8[E]=ifelse(an3<=2 && (CVM_PvalTn> 0.05),1,0)     

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
	  next	     	   
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
HN_rep=c(length(ut71[,1]),((length(ut71[,1])/mg)*100),colMeans(ut71[,c(-1,-7,-8)]),mean(biashn),mean(mse_hn1),"MES_Hn"=(((mean(ut71[,3])-mean(True_Pa))^2)+var(ut71[,3])),"CP_HN"=mean(ut71[,7])*100,DaicHN)

ut72=ut[which(HR_Cvm>=0.05),c(1,10,4,14,15,19,22,27)]
biashr=(ut72$HR_pa-ut88$TRUE_Pa);mse_hr1=((ut72$TNN2*ut72$HR_pa)^2)+(biashr)^2
HR_rep=c(length(ut72[,1]),((length(ut72[,1])/mg)*100),colMeans(ut72[,c(-1,-7,-8)]),mean(biashr),mean(mse_hr1),"MES_Hr"=(((mean(ut72[,3])-mean(True_Pa))^2)+var(ut72[,3])),"CP_HR"=mean(ut72[,7])*100,DaicHR)

ut73=ut[which(VHRT_Cvm>=0.05),c(1,11,5,16,17,20,23,28)]
biasvrht=(ut73$VHRT_pa-ut88$TRUE_Pa);mse_vrht1=((ut73$TNN3*ut73$VHRT_pa)^2)+(biasvrht)^2
VHRT_rep=c(length(ut73[,1]),((length(ut73[,1])/mg)*100),colMeans(ut73[,c(-1,-7,-8)]),mean(biasvrht),mean(mse_vrht1),"MES_Vhrt"=(((mean(ut73[,3])-mean(True_Pa))^2)+var(ut73[,3])),"CP_VHRT"=mean(ut73[,7])*100,DaicVRHT)

ut74=c("total sample","perc. fit","p-value CvM","Pa_hat","ave LL","ave UL","ave CL","ave Bias","ave New MSE","ave MSE","CP","D AIC")

ut88
ut75=data.frame(ut74,HN_rep,HR_rep,VHRT_rep)
ut75


write.csv(ut, "results/full_simulation_results.csv", row.names = FALSE)
write.csv(ut88, "results/summary_results.csv", row.names = FALSE)
write.csv(ut75, "results/model_comparison.csv", row.names = FALSE)
