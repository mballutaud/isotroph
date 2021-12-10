#---------------------------------------
#Fig 6
#--------------------------------------


## Data set ---------------------------------------------------------------------------
## load data set (Marin Leal et al., 2008) data retrieved from the authors
Data_CS<-data.frame(Date=c(0,60,120,180,240),IsoC_MPB=c(-15.11, -17.83, -15.66, -18.79, -20.82)
                    ,IsoC_Phyto=c(-21.91,-20.47,-22.65,-20.40,-21.31)
                    ,IsoC_Oyster=c(-18.04574474,-19.17558493,-18.85810134,-18.89322396,-18.60838601)
                    ,lambda=c(0.027,0.019,0.011,0.004,NA))

t_CS<-seq(0,max(Data_CS$Date),1)

# linear interpolation of sources signatures 
MPB_func<-approxfun(Data_CS$Date, Data_CS$IsoC_MPB, method='linear')
Phyto_func<-approxfun(Data_CS$Date, Data_CS$IsoC_Phyto, method='linear')
lambda_func<-approxfun(Data_CS$Date, Data_CS$lambda, method = 'const',yright = Data_CS$lambda[length(Data_CS$lambda)-1])

#parameters of trophic discrimination factor (TDF) per each source
TDF_MPB<-1 #TDF_A
TDF_PHYOM<-1 #TDF_B


best_p<-NULL
for (i in 1:4){ # 4 periods  between sampling dates
  # optimising the contributions - initial value of the optim algorithm is  the median (0.5)
  resf<-optim(par = c(0.5), fn = min_func, lower=0, upper=1, method='Brent',control=list(maxit=2000),
              data = list(t_opt=seq(Data_CS$Date[i],Data_CS$Date[i+1],1), Iso_0=Data_CS$IsoC_Oyster[i],obs=Data_CS$IsoC_Oyster[i+1],
                          lambda=lambda_func,Date=Data_CS$Date[c(i,i+1)],Phyto=Phyto_func,MPB=MPB_func,TDF_A=TDF_A,TDF_B=TDF_B))
  best_p<-cbind(best_p,resf$par)}

pMPB_DMM_func<-approxfun(cbind(Data_CS$Date[1:length(best_p)],as.numeric(best_p)), method='const',yright=best_p[length(best_p)])

Res<-RUN_TIMdyn2(t=t_CS,#time
                 state0=c(X=Data_CS$IsoC_Oyster[1]),#initial conditions
                 par=c(#two parameters lambda and Xinf i.e. forcing functions
                   lambda = lambda_func,
                   Xinf=approxfun(t, Xinf_func(t=t,sA=Phyto_func,sB=MPB_func,
                                               pB=pMPB_DMM_func,TDF_A=TDF_A,TDF_B=TDF_B))))

### SMM estimates
#instantaneous
pMPB_instant <-  SIMM(delta = Data_CS$IsoC_Oyster, sA = Data_CS$IsoC_MPB, sB = Data_CS$IsoC_Phyto,TDF_A=TDF_A,TDF_B=TDF_B)

#integration over time window
time_window <- 2*(log(2)/Data_CS$lambda) #time window is equal to twice the half-life
time_window[4]<-max(Data_CS$Date) #the last window is larger than the whole period of sampling, so it is replaced by the value for the whole period

pMPB_integr<-NULL
for (i in 1:4){ #number of sampling and number of windows
  delta_Phyto <- mean(Phyto_func(Data_CS$Date[2:5][i]-time_window[i]:Data_CS$Date[2:5][i]))
  delta_MPB <- mean(MPB_func(Data_CS$Date[2:5][i]-time_window[i]:Data_CS$Date[2:5][i]))
  pMPB_integr[i]<-SIMM(delta = Data_CS$IsoC_Oyster[i+1], sA = delta_MPB, sB = delta_Phyto,TDF_A=TDF_A,TDF_B=TDF_B)
}

p_MPB_integr_func<-approxfun(cbind(Data_CS$Date[2:5], pMPB_integr), method='const') #step function




## layout of FIGURE 6 ----------------------------------------------------------------------
## 5. graphical presentation
# graphical parameters
MPB_label  <- expression(paste(delta["s(MPB)"]))
PhyOm_label  <- expression(paste(delta["s(PhyOM)"]))
oyster_label  <- expression(paste(delta["Oyster"]))
lambda_label <-expression(paste(lambda,"(d"^-1,")", sep=" "))
lambda_leg <-expression(lambda)
pMPB_label <- expression(hat(p)["s(MPB)"])
SMM_instant_label <- expression(paste("SMM"["t"]))
SMM_integr_label <- expression(paste("SMM"[paste(Delta, "t")]))
mycol3<-c("#91064f", "#4f9106", "#064f91","#FE9A2E", "#d81c79","#ef9306","#068d91") #colors

tiff("dmm/output/Fig6.tiff", width = 789, height = 552.3, units = "px")#start of figure
par(mar=c(4,5,2,2))
# 3 plots to be combined in 2 row/ 2 columns and arranged by columns
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = FALSE))
###### A. graph of data #isotopic values of sources are corrected by the TDF values
plot(Data_CS$Date, Data_CS$IsoC_MPB+TDF_MPB, xlab = time_label, ylab = isotope_label, 
     las=1, cex.axis = 1.2, cex.lab=1.6,
     pch =16, cex=2, xlim= c(0,250),ylim=c(min(Data_CS$IsoC_Phyto),max(Data_CS$IsoC_MPB)+1.5), col=mycol3[1])
lines(Data_CS$Date, Data_CS$IsoC_MPB+TDF_MPB, col=mycol3[1], lwd=2, lty=2)
points(Data_CS$Date, Data_CS$IsoC_Phyto+TDF_PHYOM, xlab = "time", pch =16, cex=2, col=mycol3[2])
lines(Data_CS$Date, Data_CS$IsoC_Phyto+TDF_PHYOM, col=mycol3[2], lwd=2, lty=2)

points(Data_CS$Date, Data_CS$IsoC_Oyster, pch=16, cex=2, col=mycol3[3])
lines(Res$time,Res$X, lwd=3, col=mycol3[3]) #the best fit between DMM simu and oyster data
legend("topright", legend=c(MPB_label,PhyOm_label, oyster_label), 
       pch=c(16,16,16), col=mycol3[1:3],bty="n", inset = 0.002, cex=1.1)
mtext(paste("(", letters[1], ")", sep=""), side=3, line = -.7, adj=-0.2, cex=1.1, font = 2)

###### B. graph of etimated turnover
plot(t_CS,lambda_func(t_CS), xlab=time_label, ylab=lambda_label, 
     las=1, cex.axis = 1.2, cex.lab=1.6,
     pch=15, ylim=c(0,0.04), type="p", col=mycol3[4])
mtext(paste("(", letters[2], ")", sep=""), side=3, line = -.7, adj=-0.2, cex=1.1, font = 2)

###### C. graph of contribution

plot(pMPB_DMM_func(t_CS), ylim=c(0,1), xlim= c(0,250), xlab = time_label, ylab = pMPB_label,
     las=1, cex.axis = 1.2, cex.lab=1.6, col= mycol3[7], pch=21,cex=1.6)

for (i in 1:4){
  lines(c(Data_CS$Date[i+1]-time_window[i],Data_CS$Date[i+1]), c(pMPB_integr[i],pMPB_integr[i]), 
        col=mycol3[6], lwd=5 ) 
}
points(Data_CS$Date, pMPB_instant, col=mycol3[5], pch=16, cex=1.6)

legend("topleft", legend=c(SMM_instant_label, SMM_integr_label,"DMM"), pch=16, 
       col=mycol3[5:7], cex=1.1, bty="n", inset = 0.05)
mtext(paste("(", letters[3], ")", sep=""), side=3, line = -.7, adj=-0.2, cex=1.1, font = 2)

dev.off()#end of figure

## save the results of figure 6
result_MM_casestudy <- c(pMPB_instant[1:4], pMPB_integr, best_p) #last value of pMPB_instant is out of source polygon
p_estim_casestudy <- rbind(result_MM_casestudy[1:4],result_MM_casestudy[5:8],result_MM_casestudy[9:12])
colnames(p_estim_casestudy) <- c("t1", "t2", "t3", "t4")
rownames(p_estim_casestudy) <- c("SMM t", "SMM Dt", "DMM")
p_estim_casestudy
save(p_estim_casestudy, file="dmm/output/p_estim_casestudy.Rdata")
