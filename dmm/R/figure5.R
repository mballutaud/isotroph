#----------------------------------
# FIGURE 5 : 
# Bias estimate for the two statistic approaches 
# instantaneous SMM and integrated SMM
#----------------------------------

lambda_Index <- 2   # corresponds to constant lambda value at medium level, i.e. lambda = 0.02
Ratio <- c(0.1,seq(0.5,4,.5)) # a ratio is used to ensure regular distribution of the points along the x-axis
omega <- Ratio*lambda[lambda_Index]
SMM_t_label <- expression(paste("SMM"["t"]))
SMM_2w_label <- expression(paste("SMM"[paste(Delta, "t")]))
#tw<-c(2,0.5) # different values tested to relate time window to half life
tw=2

tiff("dmm/output/Fig5.tiff", width = 789, height = 552.3, units = "px") 
par(mar=c(4,5,2,2), mfrow=c(1,1))
plot(NA, NA, ylim = c(0,65), xlim = c(0,4.5), 
     xlab = ratio_lab, ylab = biais_lab, las=1, cex.lab=1.7, cex.axis=1.4)


legend("topright", legend=c(SMM_t_label, SMM_2w_label), 
       pch=16, col=mycol2, cex=1.4, inset = .0025, bg=NA, bty="n")


for (n in 1:100){# loop over the brownian sources
  for (j in 1:length(omega)){ #diet shifting : omega
    
    Results<-RUN_TIMdyn2(t,#time
                         state0=c(X=sourceB_brownian_100[ ,n][1]),#initial conditions
                         par=c(#two parameters lambda and Xinf i.e. forcing functions 
                           lambda = approxfun(t, lambda_func(t=t,lambda_0=lambda[lambda_Index],a=a[lambda_Index])), 
                           Xinf=approxfun(t, Xinf_func(t=t,sA=approxfun(t, sourceA_brownian_100[ ,n]),sB=approxfun(t, sourceB_brownian_100[ ,n]),
                                                       pB=approxfun(t, pB_func(t=t,freq=omega[j])),TDF_A=TDF_A,TDF_B=TDF_B))))
    
    # linear mixing model applied on signature at equilibrium
    pa_equilibr <- SIMM(delta = Results$Xinf, sA = sourceA_brownian_100[ ,n], sB = sourceB_brownian_100[ ,n],TDF_A=TDF_A,TDF_B=TDF_B)
    # linear mixing model applied on dynamic signature resulting from DMM 
    pA_instant <-  SIMM(delta = Results$X, sA = sourceA_brownian_100[ ,n], sB = sourceB_brownian_100[ ,n],TDF_A=TDF_A,TDF_B=TDF_B)
    #bias at t
    beta<- bias(pA_infer = pA_instant[1:500], pA_ref = pa_equilibr[1:500],
                omega = omega[j], lambda = lambda[lambda_Index], time_period = period)
    points(beta[1], beta[2]*100, col = mycol2[1], pch=16, cex=1.7)
    
    # on window proportional to lambda
    half_life <- log(2)/lambda[lambda_Index]
    
    for (m in 1:length(tw)){
      time_window <- tw[m]*round(half_life)
      # computing the moving average using time window
      diet_mean  <-Mov_Aver_W(t=t,R=Results$Xinf             ,time_window=time_window)
      sourA_mean <-Mov_Aver_W(t=t,R=sourceA_brownian_100[ ,n],time_window=time_window)
      sourB_mean <-Mov_Aver_W(t=t,R=sourceB_brownian_100[ ,n],time_window=time_window)
      
      # linear mixing model applied on signature at equilibrium and integrated
      pa_equilibr_integr <- SIMM(delta = diet_mean[time_window+1:length(t)], 
                                 sA = sourA_mean[time_window+1:length(t)], sB = sourB_mean[time_window+1:length(t)],
                                 TDF_A=TDF_A,TDF_B=TDF_B)
      # linear mixing model applied on dynamic signature from DMM and integrated
      pA_integr <-  SIMM(delta = Results$X[time_window+1:length(t)],
                         sA = sourA_mean[time_window+1:length(t)], sB = sourB_mean[time_window+1:length(t)],
                         TDF_A=TDF_A,TDF_B=TDF_B)
      #graph
      #bias on window
      beta_wind<- bias(pA_infer = pA_integr[1:(period-time_window)], 
                       pA_ref = pa_equilibr_integr[1:(period-time_window)],
                       omega = omega[j], lambda = lambda[lambda_Index], time_period = (period-time_window))
      points(beta_wind[1], beta_wind[2]*100, col = mycol2[m+1], pch=16, cex=1.7)
      
    }
  }
  print(n)}


dev.off()