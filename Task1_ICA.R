data = read.csv("/Users/harrisonknowles/Desktop/stat11_ica/group27.csv")
x = data[,1]
i = 1:299
#From the plot tau should be around 150
prob_tau = 1/298
plot(i,x) 

integranda <- function(sigsq,tau) {
  u<-1
  for(i in 1:tau){
    u = u * 1/(sqrt(2*pi*sigsq)) * exp(-1/2*(x[i]/sqrt(sigsq))^2)
  }
  return(u)
}

integranda2 <- function(sigsq,tau) {
  u<-1
  for(i in (tau:298)){
    u = u * 1/(sqrt(2*pi*sigsq)) * exp(-1/2*(x[i]/sqrt(sigsq))^2)
  }
  return(u)
}

integrandb <- function(sigsq){
  v = (1/(sigsq))^2 * exp(-1/sigsq)
  return(v)
}

integrand1 <- function(sigsq,tau){
  int = integranda(sigsq,tau)*integrandb(sigsq)
  return(int)
}

integrand2 <- function(sigsq,tau){
  int = integranda2(sigsq,tau)*integrandb(sigsq)
  return(int)
}

probs<- rep(0,298)
j<-1
while(j < length(probs)){
  probs[j]=(integrate(integrand1,0,Inf,tau=j)$value)*(integrate(integrand2,0,Inf,tau=j)$value)
  j = j+1
}
probs = (probs*prob_tau) / sum(probs)
sum(probs)
max(probs)
which.max(probs)
probs
plot(i,x)
abline(v=145, col="blue")

