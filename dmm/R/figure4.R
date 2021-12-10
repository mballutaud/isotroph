#----------------------------------
# FIGURE 4 : 
# Estimated contributions of source A to a consumer's diet compared to reference.
# Contribution are inferred from isotopic composition of consumer simulated 
# using DMM and forcing food sources over time.
#----------------------------------
tiff("dmm/output/Fig4.tiff", width = 789, height = 552.3, units = "px") 
lambda_Index <- 2   # corresponds to constant lambda value at medium level, i.e. lambda = 0.02
SMM_t_label <- expression(paste("SMM"["t"]))
SMM_2w_label <- expression(paste("SMM"[paste(Delta, "t")]))
# parameters for integrated SMM over a time window
half_life <- log(2)/lambda[lambda_Index]
time_window <- round(2*half_life)


par(mfrow=c(2,3), mar=c(4,6,3,1)) # multi-panels graph

### FIG 4a,b,c : Comparison between instantaneous SMM and DMM 
# FIG 4 a,b : Food sources constant
for (j in 1:2){ # omega for 1 and 4 diet-switches
  
  Results<-RUN_TIMdyn2(t,                       #time
                       state0=c(X=sB_t[1]),     #initial conditions
                       par=c(                   #two parameters lambda and Xinf i.e. forcing functions 
                         lambda = approxfun(t, lambda_func(t=t,lambda_0=lambda[lambda_Index],a=a[lambda_Index])), 
                         Xinf=approxfun(t, Xinf_func(t=t,sA=approxfun(t, sA_t),sB=approxfun(t, sB_t),
                                                     pB=approxfun(t, pB_func(t=t,freq=omega[j])),TDF_A=TDF_A,TDF_B=TDF_B))))
  
  pA_instant <-  SIMM(delta = Results$X, sA = sA_t, sB = sB_t,TDF_A=TDF_A,TDF_B=TDF_B) # linear mixing model applied on dynamic signature resulting from DMM
  pA_equilibr <- SIMM(delta = Results$Xinf, sA = sA_t, sB = sB_t,TDF_A=TDF_A,TDF_B=TDF_B) # linear mixing model applied on the signature at equilibrium
  
  plot(t, pA_instant, col=mycol2[1], lwd=3, lty=1, type = "l", las=1, cex.axis = 1.4, cex.lab =1.4,
       xlab = "", ylab = "", xlim = c(0,500), ylim = c(0,1))
  lines(t, pA_equilibr, col=mycol2[3], lwd=3, lty=2) 
  mtext(paste("(", letters[j], ")", sep=""), side=3, line = -1, adj=-0.25, cex=1.1, font = 2)
  if (j==1){
    legend("bottomright", legend = SMM_t_label, lty = 1, col= mycol2[1], lwd = 2, bty = "n", cex=1.4, inset = 0.05)
  }
}

# FIG 4c : Food sources Brownian 
omega_index = 2 # for browninan trajectories, only 4 diet switches are represented

Results<-RUN_TIMdyn2(t,#time
                     state0=c(X=sB_t_BR[1]),#initial conditions
                     par=c(#two parameters lambda and Xinf i.e. forcing functions 
                       lambda = approxfun(t, lambda_func(t=t,lambda_0=lambda[lambda_Index],a=a[lambda_Index])), 
                       Xinf=approxfun(t, Xinf_func(t=t,sA=approxfun(t, sA_t_BR),sB=approxfun(t, sB_t_BR),
                                                   pB=approxfun(t, pB_func(t=t,freq=omega[omega_index])),TDF_A=TDF_A,TDF_B=TDF_B))))

# linear mixing model applied on dynamic signature resulting from DMM
pA_instant <-  SIMM(delta = Results$X, sA = sA_t_BR, sB = sB_t_BR, TDF_A=TDF_A,TDF_B=TDF_B)
# linear mixing model applied on the signature at equilibrium
pA_equilibr <- SIMM(delta = Results$Xinf, sA = sA_t_BR, sB = sB_t_BR,TDF_A=TDF_A,TDF_B=TDF_B)

plot(t, pA_instant, col=mycol2[1], lwd=3, lty=1, type = "l", las=1, cex.axis = 1.4, cex.lab =1.4,
     xlab = "", ylab = "", xlim = c(0,500), ylim = c(0,1))
lines(t,pA_equilibr, col=mycol2[3], lwd=3, lty=2) 
mtext(paste("(", letters[3], ")", sep=""), side=3, line = -1, adj=-0.25, cex=1.1, font = 2)

### FIG 4d,e,f : Comparison between integrated SMM and DMM 
# FIG 4 d,e : Food sources constant
for (j in 1:2){ # omega for 1 and 4 diet-switches
  
  Results<-RUN_TIMdyn2(t,#time
                       state0=c(X=sB_t[1]),#initial conditions
                       par=c(#two parameters lambda and Xinf i.e. forcing functions 
                         lambda = approxfun(t, lambda_func(t=t,lambda_0=lambda[lambda_Index],a=a[lambda_Index])), 
                         Xinf=approxfun(t, Xinf_func(t=t,sA=approxfun(t, sA_t),sB=approxfun(t, sB_t),
                                                     pB=approxfun(t, pB_func(t=t,freq=omega[j])),TDF_A=TDF_A,TDF_B=TDF_B))))
  
  # perform moving average on forcing diet and food sources signatures
  diet_integr  <-Mov_Aver_W(t=t,R=Results$Xinf,time_window=time_window)
  sourA_integr <-Mov_Aver_W(t=t,R=sA_t    ,time_window=time_window)
  sourB_integr <-Mov_Aver_W(t=t,R=sB_t    ,time_window=time_window)
  
  # linear mixing model applied on the signature at equilibrium integrated over time window:
  pA_equilibr_integr <- SIMM(delta = diet_integr[time_window+1:length(t)], 
                             sA = sourA_integr[time_window+1:length(t)], sB = sourB_integr[time_window+1:length(t)],
                             TDF_A=TDF_A,TDF_B=TDF_B)
  # linear mixing model applied on dynamic signature resulting from DMM and integrated over time window:
  pA_integr <-  SIMM(delta = Results$X[time_window+1:length(t)],
                     sA = sourA_integr[time_window+1:length(t)], sB = sourB_integr[time_window+1:length(t)],
                     TDF_A=TDF_A,TDF_B=TDF_B)
  
  plot(t[time_window+1:length(t)],pA_integr, col=mycol2[2], lwd=3, lty=1, type = "l", las=1, cex.axis = 1.4, cex.lab =1.4,
       xlab = "", ylab = "", xlim = c(0,500), ylim = c(0,1))
  if (j==1){
    legend("bottomleft", lty = c(2,1), legend = c(pA_label,pA_res), bty = "n", cex=1.2)
    legend("bottomright", legend = SMM_2w_label, lty = 1, col= mycol2[2], lwd = 2, bty = "n", cex=1.4, inset = 0.025)
  }
  lines(t[time_window+1:length(t)], pA_equilibr_integr,col=mycol2[3], lwd=3, lty=2)
  mtext(paste("(", letters[3+j], ")", sep=""), side=3, line = -1, adj=-0.25, cex=1.1, font = 2)
}  


# FIG 4f : Food sources Brownian 
omega_index = 2 # for browninan trajectories, only 4 diet switches are represented

Results<-RUN_TIMdyn2(t,#time
                     state0=c(X=sB_t_BR[1]),#initial conditions
                     par=c(#two parameters lambda and Xinf i.e. forcing functions 
                       lambda = approxfun(t, lambda_func(t=t,lambda_0=lambda[lambda_Index],a=a[lambda_Index])), 
                       Xinf=approxfun(t, Xinf_func(t=t,sA=approxfun(t, sA_t_BR),sB=approxfun(t, sB_t_BR),
                                                   pB=approxfun(t, pB_func(t=t,freq=omega[omega_index])),TDF_A=TDF_A,TDF_B=TDF_B))))


# perform moving average on forcing diet and food sources signatures   
diet_integr  <-Mov_Aver_W(t=t,R=Results$Xinf,time_window=time_window)
sourA_integr <-Mov_Aver_W(t=t,R=sA_t_BR ,time_window=time_window)
sourB_integr <-Mov_Aver_W(t=t,R=sB_t_BR ,time_window=time_window)

# linear mixing model applied on the signature at equilibrium integrated over time window:
pA_equilibr_integr <- SIMM(delta = diet_integr[time_window+1:length(t)], 
                           sA = sourA_integr[time_window+1:length(t)], sB = sourB_integr[time_window+1:length(t)],
                           TDF_A=TDF_A,TDF_B=TDF_B)

# linear mixing model applied on dynamic signature resulting from DMM and integrated over time window:
pA_integr <-  SIMM(delta = Results$X[time_window+1:length(t)],
                   sA = sourA_integr[time_window+1:length(t)], sB = sourB_integr[time_window+1:length(t)],
                   TDF_A=TDF_A,TDF_B=TDF_B)

plot(t[time_window+1:length(t)], pA_integr, col=mycol2[2], lwd=3, lty=1, type = "l", las=1, cex.axis = 1.4, cex.lab =1.4,
     xlab = "", ylab = "", xlim = c(0,500), ylim = c(0,1))

lines(t[time_window+1:length(t)], pA_equilibr_integr,col=mycol2[3], lwd=3, lty=2)

# adding axis labels and legend
mtext(pA_label, side=2, line=-2, cex=1.5, outer = TRUE)
mtext(time_label, side=1, line=-1, cex=1.5, outer = TRUE)

mtext(paste("(", letters[6], ")", sep=""), side=3, line = -1, adj=-0.25, cex=1.1, font = 2)

omega_one <- expression(paste(omega["one switch"]))
omega_four <- expression(paste(omega["four switches"]))
delta_const <- expression(paste(delta["constant"]))
delta_brown<- expression(paste(delta["Brownian"]))

mtext(omega_one, side=3, line=-2, cex=1.2, adj = 0.12, outer = TRUE)
mtext(delta_const, side=3, line=-2, cex=1.2, adj = 0.22, outer = TRUE)

mtext(omega_four, side=3, line=-2, cex=1.2, adj = 0.5, outer = TRUE)
mtext(delta_const, side=3, line=-2, cex=1.2, adj = 0.6, outer = TRUE)

mtext(omega_four, side=3, line=-2, cex=1.2, adj = 0.85, outer = TRUE)
mtext(delta_brown, side=3, line=-2, cex=1.2, adj = 0.95, outer = TRUE)

dev.off()