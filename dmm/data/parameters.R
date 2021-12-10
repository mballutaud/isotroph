####################################################
####     GLOBAL PARAMETERS    ----------------------
####################################################

period <-500 #time period of simulation (d)
t <- seq(0,period,1) # output times for the ODE (d)

nswitch <- c(1,4,8) # number of switches between food sources
omega <- 1/period*nswitch # frequency of diet shift

# modalities for lambda
# using the generic function lambda_func, lambda can be either set constant (with a = 0) or dynamic (with a>0)
lambda <- c(0.002, 0.02, 0.2, 0.2)
a<-c(0,0,0, 0.01)
# here, 4 modalities : constant low, constant medium, constant high, and dynamicly decreasing from 0.2


# FOOD SOURCES -----------------------------------------------------------------------------------------
# isotope signatures of food sources A and B 
# could be constant (sA_t and sB_t) or variable following a brownian trajectory (sA_t_BR and sB_t_BR)

#load("random-sources-100simu.Rdata") # matrix of Brownian trajectories - see script product-brownian-sources.R to produce them
brownian_index = 2 # a number between 1 and 100 to select one trajectory among the 100 available - published figure was realized with IDbrowninan = 2

# source A
# if constant
sA_cst <- 0 # isotope ratio value 
sA_t <- rep(sA_cst, length(t))
# if brownian
sA_t_BR <- sourceA_brownian_100[,brownian_index] 

# source B
# if constant
sB_cst <- 10 # isotope ratio value 
sB_t <- rep (sB_cst, length(t))
# if brownian
sB_t_BR <- sourceB_brownian_100[,brownian_index] 

#Trophic Discrimination Factors, possibly different across food sources A and B
TDF_A<-1
TDF_B<-1



# FIGURES: COMMON PARAMETERS AND LABELS -----------------------------------------------------------

# creating labels using greek alphabet for delta and lambda
# - for axes
isotope_label <- expression(paste(delta," (\u2030)", sep=" "))
time_label <- paste("Time", "(d)")
# - for legend
L_lambda <- expression(paste(delta["c"]," (",lambda, " =  2",".","10"^"-3"," d"^"-1",")"))
M_lambda <- expression(paste(delta["c"]," (",lambda, " =  2",".","10"^"-2"," d"^"-1",")"))
H_lambda <- expression(paste(delta["c"]," (",lambda, " =  2",".","10"^"-1"," d"^"-1",")"))
onto_lambda <- expression(paste(delta["c"]," (","ontogenetic ", lambda,")"))
diet_label <- expression(paste(delta["d"], sep=" "))
sour_label <- expression(paste(delta["s(i)"], sep=" "))
ratio_lab <- expression(paste(omega, "/", lambda))
biais_lab <- expression(paste(beta, "(%)"))
pA_label <- expression(paste("p"["s(a)"]))
#pA_res <- expression(paste("p'"["s(a)"]))
pA_res <- expression(hat(p)["s(a)"])

mycol<-c("#FF3366", "#3399CC", "#66CC33", "#FE9A2E","#333333","#848484")
mycol2<- c("#d81c79","#ef9306","#068d91","#78909c") #ffab40 orange #0097a7 aqua
