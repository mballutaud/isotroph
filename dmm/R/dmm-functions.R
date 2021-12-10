###### DMM and associated functions
# previous function by Marine
#pA_func <- function(t,freq) {round((cos(pi*freq*(t+1/(2*freq)))+1)/2)} #binary value as periodic square wave

# (Equations 7 in Ballutaud et al.)
# gives values of 0 and 1 alternatively
pB_func <- function(t,freq) {1-ceiling(sin(pi*freq*t)/2)} 

# IsoDietfunc computes the isotopic signature of the diet composed of 2 food sources A and B  
IsoDiet_func <- function(t,sA,sB,pB){(1-pB)*sA+pB*sB}

# Function providing signature of consumer at equilibrium 
Xinf_func <- function(t,sA,sB,pB,TDF_A,TDF_B){(1-pB(t))*(sA(t)+TDF_A)+(pB(t))*(sB(t)+TDF_B)}


# Function implementing the linear mixing-model for 2 sources (Phillips & Gregg 2003)
# with delta the signature of the consumer, Sa and sB the signatures of the two sources
# return proportion of source A
SIMM <- function(delta, sA, sB, TDF_A, TDF_B){ 
  pA <- (delta-(sB+TDF_B))/(sA+TDF_A-(sB+TDF_B))
  return(pA)
}


#numerical solution of ODE
#with state0 the initial conditions 
#and par composed of lambda and Xinf forcing functions 
RUN_TIMdyn2<-function(t,state0,par){ 
  out<-as.data.frame(lsoda(state0,t,ODE_TIMdyn2,par)) # lsoda is an integrator to solve ODE
  return(out)}

# differential equation of the dynamic mixing model
ODE_TIMdyn2<-function(t,state,parameters){ #time incorporation model : ODE
  with(as.list(c(state,parameters)),{  # unpack the state variables, parameters
    #state variable X(t) : the response
    #Forcing variables lambda(t) and Xinf(t)
    dX <- lambda(t)*(Xinf(t)-X)
    # the outputs
    list(c(dX = dX), c(Xinf = Xinf(t), lambda = lambda(t))) #the response is saved
  })
}

# function computing the difference between DMM and observed data
min_func<-function(data,par){
  Res<-RUN_TIMdyn2(t=data$t_opt,#time
                   state0=c(X=data$Iso_0),#initial conditions
                   par=c(#two parameters lambda and Xinf i.e. forcing functions 
                     lambda = data$lambda, 
                     Xinf=approxfun(t, Xinf_func(t=t,sA=data$Phyto,sB=data$MPB,
                                                 pB=approxfun(cbind(data$Date,c(par,par)), method='const',yright = par),TDF_A=data$TDF_A,TDF_B=data$TDF_B))))
  cost<-abs(data$obs-Res$X[length(Res$X)])
  return(cost)}


# Function for getting lambda values for all time steps (vector t), either constant at lambda_0 (if a=0) or decreasing with time from lambda_0 (with a>0)
# (Equation 9 in Ballutaud et al.) 
lambda_func <- function(t,lambda_0,a){ 
  lambda <- lambda_0*exp(-a*t)
  return(lambda)
}


# function computing moving average of R over a specified time window
Mov_Aver_W <- function(t,R,time_window){
  MAW<-NULL
  for (s in 1:length(t)){MAW[s+time_window] <- mean(R[s:(s+time_window)])}
  return(MAW)}

# function computing the bias between pa_infer and pA_reference
bias <- function(pA_infer, pA_ref, time_period, omega, lambda){ #bias estimates
  ratio <- round(omega/lambda, digits = 2)
  B <- round(sum(abs(pA_infer-pA_ref))/time_period, digits = 3)
  return(c(ratio,B))
}