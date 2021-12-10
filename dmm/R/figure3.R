#----------------------------------
# FIGURE 3 : 
# Simulated isotopic values of consumer (delta_c(t)) over 500 d, for
# contrasted isotopic turnover rates (constant and ontogenetic) and experiencing a
# variable diet (delta_d(t) dashed line) resulting from four diet-switches
# between two food sources as simulated with Brownian trajectories
#----------------------------------

tiff("dmm/output/Fig3.tiff", width = 789, height = 552.3, units = "px") 
par(mar=c(4,4,2,2),mfrow=c(1,1))  

plot(NULL,  ylab = "", xlab = "", xlim=c(0,500),ylim = c(-5,15), las=1, cex.axis = 1.2) 

# compute the isotope signature of the diet
# omega[2] corresponds to 4 diet switches
Diet<-IsoDiet_func(t=t,sA=sA_t_BR,sB=sB_t_BR,pB=pB_func(t=t,freq=omega[2]))+1 #+1 is the TDF

# loop over the 4 modalities for lambda : constant low, constant medium, constant high, and dynamicaly decreasing from 0.002
for (k in 1:4){ #turnover
  
  Results<-RUN_TIMdyn2(t,#time
                       state0=c(X=sB_t[1]+TDF_B),    #initial conditions
                       par=c(    #two parameters lambda and Xinf i.e. forcing functions 
                         lambda = approxfun(t, lambda_func(t=t,lambda_0=lambda[k],a=a[k])), 
                         Xinf=approxfun(t,Xinf_func(t=t,sA=approxfun(t, sA_t_BR),sB=approxfun(t, sB_t_BR),
                                                    pB=approxfun(t, pB_func(t=t,freq=omega[2])),TDF_A=TDF_A,TDF_B=TDF_B))))
  
  lines(Results$t, Results$X,lty=1, lwd=3, col= mycol[k])
}

# adding the line of the diet signatures
lines(t,Diet, lwd=2, lty=2, col=mycol[5])
# adding the lines of the two food sources signatures, identical for all experiments
lines(t, sA_t_BR, lwd=2, lty=3, col=mycol[6])
lines(t, sB_t_BR, lwd=2, lty=3, col=mycol[6])
# adding axis legend
mtext(isotope_label, side=2, line=2, cex=1.6)
mtext(time_label, side=1, line=2.5, cex=1.6)

legend("bottomright", legend=c(L_lambda, M_lambda, H_lambda, onto_lambda, diet_label, sour_label), col=mycol,
       lty = c(1,1,1,1,2,3), lwd=2, inset = 0.03, cex=1.1, bty = "n")#,y.intersp=0.7)

dev.off()