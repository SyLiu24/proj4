
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}
rb1<-function(th){
  th^3
}
gb1<-function(th){
  2*th^2
}
hb1<-function(th){
  4*th
}
rb2<-function(th){
  log(th,base=exp(1))
}
gb2<-function(th){
  1/th
}
hb2<-function(th){
  -1/th^2
}
newt(c(0,0),rb,gb,hb,k=2)




nll<-function(th,t,y){
  mu<-th[1]*exp(th[2]*t)
  -sum(dpois(y,mu,log=TRUE))
}
# gll<-function(th,t,y){
#   c(sum(-y/th[1]+exp(t*th[2])),sum(-y*t+t*th[1]*exp(t*th[2])))
# }
gll <- function(theta,t,y) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t) ## avoid computing twice
  -c(sum(y)/alpha - sum(ebt), ## -dl/dalpha
     sum(y*t) - alpha*sum(t*ebt)) ## -dl/dbeta
} ## gll

hll <- function(theta,t,y) {
  ## Hessian of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t) ## avoid computing twice
  H <- matrix(0,2,2) ## matrix for Hessian of -ve ll
  H[1,1] <- sum(y)/alpha^2
  H[2,2] <- alpha*sum(t^2*ebt)
  H[1,2] <- H[2,1] <- sum(t*ebt)
  H
} ## hll




t=1:13
y=c(12,14,33,50,67,74,123,141,165,204,253,246,240)
fit <- optim(c(10,.1),nll,gr=gll,y=y,t=t,method="BFGS",hessian=TRUE)
fit