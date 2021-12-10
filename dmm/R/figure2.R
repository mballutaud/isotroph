#----------------------------------
# FIGURE 2 : 
# Simulated isotopic values of consumer (delta_c(t)) during 500 d, for
# contrasted isotopic turnover rates (Lambda, solid lines) and experiencing a
# variable diet (delta_d(t) dashed line) resulting from one or four diet-switches
# between two food sources
#----------------------------------
tiff("dmm/output/Fig2.tiff", width = 789, height = 552.3, units = "px") #Minimum width in pixels for plos one = 789 and h =789*7/10
par(mar=c(4,4,3,1), mfrow=c(1,2))

leg_fig2<-c("(a) One diet switch","(b) Four diet switches") # legend for the two panels

# Generating results and plot simultaneously for 2 factors (Frequency of diet shift and Isotopic turnover rate), one loop each
for (j in 1:2){ #diet shift
  plot(NULL, ylab = "", xlab = "", xlim=c(0,period), ylim = c(min(sA_cst, sB_cst),max(sA_cst, sB_cst)+max(TDF_A,TDF_B)), las=1, cex.axis = 1.2)
  
  # compute the isotope signature of the diet
  Diet<-IsoDiet_func(t=t,sA=sA_t,sB=sB_t,pB=pB_func(t=t,freq=omega[j]))+1 #plus the TDF value
  
  for (k in 1:3){ #Isotopic turnover rate
    
    # run the DMM
    Results<-RUN_TIMdyn2(t,  #time
                         state0=c(X=sB_t[1]+TDF_B),   #initial conditions = ISOTOPIC DIET
                         par=c(    #two parameters lambda and Xinf i.e. forcing functions 
                           lambda = approxfun(t, lambda_func(t=t,lambda_0=lambda[k],a=a[k])), 
                           Xinf=approxfun(t,Xinf_func(t=t,sA=approxfun(t, sA_t),sB=approxfun(t, sB_t),
                                                      pB=approxfun(t, pB_func(t=t,freq=omega[j])),TDF_A=TDF_A,TDF_B=TDF_A))))
    
    lines(Results$t, Results$X, lty=1, lwd=3, col= mycol[k])
  }
  
  # adding the line of the diet signatures
  lines(t,Diet, lwd=2, lty=2, col=mycol[5])
  # adding axis legend and letter of the panel
  mtext(leg_fig2[j], side=3, line =.75, adj=-0.1, cex=1.4, font = 2)
  mtext(isotope_label, side=2, line=2, cex=1.6)
  mtext(time_label, side=1, line=2.5, cex=1.6)
  if (j==1){
    legend("topright",   # Coordinates (x also accepts keywords)
           legend=c(L_lambda, M_lambda, H_lambda,diet_label), # Vector with the name of each group
           col=c(mycol[1:3], mycol[5]), # Color of lines or symbols
           lty = c(1,1,1,2), lwd=2,         # Line type and width
           bty = "n",        # Box type (bty = "n" removes the box)
           cex = 1.1,          # Legend size
           inset = 0.01)
  }
}
dev.off()