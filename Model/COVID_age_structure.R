#This script contains an age-structure model with heterogeneous mixing
library(deSolve)
library(tidyverse)

##############################
#Read in pcms - I'll deliver these in better formats

c_in <- read_csv("Data/all_pcms.csv") 

C <- add_column(c_in,age=rep(names(c_in)[1:3],length.out=nrow(c_in))) %>%
  group_by(age) %>%
  summarize_at(vars(1:3),sum) %>%
  select(-1) %>%
  as.matrix()


#############################


seir_age <- function(t, x, vparameters){
  ncompartment = 4
  nage = length(x)/ncompartment
  S    = as.matrix(x[1:nage])
  E    = as.matrix(x[(nage+1):(2*nage)])
  I    = as.matrix(x[(2*nage+1):(3*nage)])
  R    = as.matrix(x[(3*nage+1):(4*nage)])
  
  I[I<0] = 0
  with(as.list(vparameters),{
    # note that because S, I and R are all vectors of length nage, so will N,
    # and dS, dI, and dR
    N = S+E+I+R
    dS = -beta * S * C %*% I/N
    dE = beta * S * C %*% I/N - E/alpha
    dI = E/alpha - gamma*I
    dR = gamma*I
    # remember that you have to have the output in the same order as the model
    # compartments are at the beginning of the function
    out=c(dS,dE,dI,dR)
    list(out)
  })
}

npop = 10000000
f = c(0.3,0.45,0.25) # two age classes, with 25% kids, and 75% adults
N =  npop*f      # number in each age class

nage = length(f)
I_0    = rep(1,nage) # put one infected person each in the kid and adult classes
S_0    = N-I_0
E_0    = rep(0,nage)
R_0    = rep(0,nage)

gamma = 1/5        # recovery period of influenza in days^{-1}
R0    = 2.68        # R0 of covid
alpha = 5.1 ## incubation period





lcalculate_transmission_probability = 1 # if this is 1, then calculate the transmission probability from R0

if (lcalculate_transmission_probability==1){
  eig = eigen(C)
  # reverse engineer beta from the R0 and gamma 
  beta = R0*gamma/max(Re(eig$values))  
  beta = beta
}else{
  beta = c(0.05)
}
vparameters = c(gamma=gamma,beta=beta,C=C)
inits = c(S=S_0,E=E_0,I=I_0,R=R_0)

##################################################################################
# let's determine the values of S,I and R at times in vt
##################################################################################
vt = seq(0,200,1)  
mymodel_results = as_tibble(lsoda(inits, vt, seir_age, vparameters))


mymodel_results %>%
  mutate_all(as.numeric) %>%
  select(time,contains("I")) %>%
  pivot_longer(-time,names_to = "class",values_to = "values") %>%
  ggplot(aes(x=time,y=values,color=class)) +
  geom_line()

