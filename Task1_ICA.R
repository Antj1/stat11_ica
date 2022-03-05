data = read.csv("/Users/harrisonknowles/Desktop/stat11_ica/group27.csv")
x = data[,1]
i = 1:299
#From the plot tau should be around 150
prob_tau = 1/299
plot(i,x) 

integranda <- function(sig,tau) {
  for(i in range(tau)){
    u = 1
    u = u * 1/(sqrt(2*pi)*sig) * exp(-1/2*(x[i]/sig)^2)
  }
  return(u)
}

integranda2 <- function(sig,tau) {
  for(i in range(tau,299)){
    u = 1
    u = u * 1/(sqrt(2*pi)*sig) * exp(-1/2*(x[i]/sig)^2)
  }
  return(u)
}

integrandb <- function(sig){
  v = (1/(sig^2))^2 * exp(-1/(sig)^2)
  return(v)
}

integrand1 <- function(sig,tau){
  int = integranda(sig,tau)*integrandb(sig)
  return(int)
}

integrand2 <- function(sig,tau){
  int = integranda2(sig,tau)*integrandb(sig)
  return(int)
}

probs<- rep(1,299)
j<-1
while(j < 300){
  probs[j]=(integrate(integrand1,0,Inf,tau=j)$value)*(integrate(integrand2,0,Inf,tau=j)$value)
  j = j+1
}
probs = probs / sum(probs)
sum(probs)
