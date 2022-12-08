
#####PREDATOR PREY SIMULATIONS######

##M DeSiervo mdesierv@uwyo.edu###

## Dec 2022 ###

#In this exercise, we will model and simulate data for the classic Lotka-Volterra predator-prey
#model. We will build intuition by first developing a model for prey with density-dependent intraspecific growth## 
## we will then build the classic Lotka volterra predator prey model, and add a bit more complexity by adding back in 
##intraspecific density dependence in prey population, and a type 2 functional response for predator###


#install.packages(deSolve) 
library(deSolve) 

############ PREY MODEL WITH INTRASPECIFIC DENSITY DEPENDENCE######

## Prey = N ### 
# ----------------------------------------------------------------------------------------------------
# Define Function to run model with prey species for overall dynamics
preymodel <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN1 = r*N*(1-N/K)  ## prey w/ density dep pop growth###
    return(list(dN1))
  })
}
# ----------------------------------------------------------------------------------------------------


##Define parameters##

###r growth rate of prey##

r=0.2

## K carrying capacity of prey##

K= 100


#########

State <- c(N=1) # starting conditions for simulation## ##

pars <- c(r, K) ##parameters for simulation ##

end_time <- 200  ## run the model for 200 times steps##
times <- seq(1, end_time, by = 1) 

######

# solve system of equations using the ode method with lsoda
solnprey1=lsoda(y=State, times=times, func=preymodel , parms=pars, hmax=0.1)

###Plot the output###
plot(solnprey1[,2], type="l", col="black", lwd=2, xlab="Time", ylab="Number of Individuals", xlim=c(0,200), ylim=c(0,100))
lines(solnprey1[,2], lwd=2, col="red")
legend(x=80, y=90, c("N = Prey"), col=c("red"), lwd=2)

##### how does changing estimates for r and K change this model?######



############ PREDATOR - PREY MODEL ###### Classic Lotka Volterra model####


# ----------------------------------------------------------------------------------------------------
# Function to run model with 2 species for overall dynamics
Predpreymodel1 <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    ##prey ##
    dN = r*N-a*N*P    
    
    ##predator##
    dP = b*a*N*P- d*P
    
    
    return(list(c(dN, dP)))
  })
}
# ----------------------------------------------------------------------------------------------------


# parameters


###r growth rate of prey##

r=0.2

### a attack rate of predator on prey##

a = 0.05

###b conversion efficiency prey into predator babies##

b = 0.2

## d = predator death rate

d = 0.2


##

State <- c(N=50, P=2) # starting conditions this time start with 50 prey and 2 predator##

pars <- c(r, a, b, d) ##parameters##

end_time <- 200
times <- seq(1, end_time, by = 1) 

#

# solve system of equations using the ode method with lsoda

solnpredprey1=lsoda(y=State, times=times, func=Predpreymodel1 , parms=pars, hmax=0.1)

### plot the output for predator prey model###

plot(solnpredprey1[,2], type="l", col="black", lwd=2, xlab="Time", ylab="Number of Individuals", xlim=c(0,200), ylim=c(0,60))
lines(solnpredprey1[,2], lwd=2, col="red")
lines(solnpredprey1[,3], lwd=2, col="blue")
legend(x=160, y=45, c("N = Prey", "P = Predator"), col=c("red", "blue"), lwd=2)

###populations cycle around unstable equiblrium! ###

####what happens if we increase conversion efficiency? ###



############ PREDATOR - PREY MODEL w/ prey DD ####


# ----------------------------------------------------------------------------------------------------
# Function to run model with 2 species for overall dynamics
Predpreymodeldd <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    ##prey ##
    dN =  r*N*(1-N/K)-a*N*P     
    
    ##predator##
    dP = b*a*N*P- d*P
    
    
    return(list(c(dN, dP)))
  })
}
# ----------------------------------------------------------------------------------------------------


# parameters


###r growth rate of prey##

r=0.2

### a attack rate of predator on prey##

a = 0.05

###b conversion efficiency prey into predator babies##

b = 0.2

## d = predator death rate

d = 0.2


##

State <- c(N=50, P=2) # starting conditions this time start with 50 prey and 2 predator##

pars <- c(r, a, b, d) ##parameters##

end_time <- 200
times <- seq(1, end_time, by = 1) 

#

# solve system of equations using the ode method with lsoda

solnpredprey2=lsoda(y=State, times=times, func=Predpreymodeldd , parms=pars, hmax=0.1)

### plot the output for predator prey model###

plot(solnpredprey2[,2], type="l", col="black", lwd=2, xlab="Time", ylab="Number of Individuals", xlim=c(0,200), ylim=c(0,60))
lines(solnpredprey2[,2], lwd=2, col="red")
lines(solnpredprey2[,3], lwd=2, col="blue")
legend(x=160, y=45, c("N = Prey", "P = Predator"), col=c("red", "blue"), lwd=2)


###Damped oscillations until reaching stable equilibrium###



############ PREDATOR - PREY MODEL w/ prey DD and Type II functional response ####


# ----------------------------------------------------------------------------------------------------
# Function to run model with 2 species for overall dynamics
Predpreymodelddtype2 <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    ##prey ##
    dN =  r*N*(1-N/K)-a*N*P     
    
    ##predator##
    dP = b*(a*N/(1+a*h*N))- d*P ## add in type 2 functional response### 
    
    
    return(list(c(dN, dP)))
  })
}
# ----------------------------------------------------------------------------------------------------


# parameters


###r growth rate of prey##

r=0.3

### a attack rate of predator on prey##

a = 0.05

###b conversion efficiency prey into predator babies##

b = 0.3

###h = handling time###

h =0.1

## d = predator death rate

d = 0.1


##

State <- c(N=50, P=2) # starting conditions this time start with 50 prey and 2 predator##

pars <- c(r, a, b, d) ##parameters##

end_time <- 200
times <- seq(1, end_time, by = 1) 

#

# solve system of equations using the ode method with lsoda

solnpredprey3=lsoda(y=State, times=times, func=Predpreymodelddtype2 , parms=pars, hmax=0.1)

### plot the output for predator prey model###

plot(solnpredprey3[,2], type="l", col="black", lwd=2, xlab="Time", ylab="Number of Individuals", xlim=c(0,200), ylim=c(0,60))
lines(solnpredprey3[,2], lwd=2, col="red")
lines(solnpredprey3[,3], lwd=2, col="blue")
legend(x=160, y=45, c("N = Prey", "P = Predator"), col=c("red", "blue"), lwd=2)


##Stable equilibrium at a much lower abundance of prey###
