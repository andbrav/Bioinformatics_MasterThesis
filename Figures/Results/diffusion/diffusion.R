library(ReacTran)
library(ggplot2)
require(deSolve)

N=45000000
A=505990
D=N/A
beta=0.45
k= 150000/A
rmu=0.021
n=10
b=0.03
D1=10
D2=5



Gsize<- 1000
Grid <- setup.grid.1D(x.up = 0, x.down = A, N=Gsize)

Sini<-rep((N-1)/A,times=Gsize)
Iini<-rep(0,times=Gsize)
Iini[Gsize%/%2]<- 1/(A/Gsize)
Rini<-rep(0,times=Gsize)
initial_values<-c(Sini,Iini,Rini)

q <- function(X,b2,n2,k2) {b2*k2**n2/(k2**n2+ X**n2)}
SCIR1D <- function(t,state,params){
  S <- state[1:Gsize]
  I <- state[(Gsize+1):(2*Gsize)]
  
  dS <- -b*S + q(I,b,n,k)*D + tran.1D(C=S,C.up=S[1],C.down=S[Gsize],D=D1,dx=Grid,AFDW=0.5)$dC
  dI <- (beta*S*I)/D - rmu*I + tran.1D(C=I,D=D2,dx=Grid,AFDW=0.5)$dC
  dR <- rmu*I
  
  list(c(dS,dI,dR))
}



times <- seq(1, 10000, by=1/10)
parameters <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu,D=D,D1=2,D2=1)
print(system.time(
out1 <- ode.1D(y = initial_values, func = SCIR1D,
              times = times, parms = parameters, nspec = 2,
              names = c("S", "I"), dimens = Gsize, method="euler")
))

times <- seq(1, 10000, by=1/10)
parameters <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu,D=D,D1=2,D2=1)
print(system.time(
  out2 <- ode.1D(y = initial_values, func = SCIR1D,
                 times = times, parms = parameters, nspec = 2,
                 names = c("S", "I"), dimens = Gsize, method="euler")
))

times <- seq(1, 10000, by=1/10)
parameters <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu,D=D,D1=2,D2=1)
print(system.time(
  out3 <- ode.1D(y = initial_values, func = SCIR1D,
                 times = times, parms = parameters, nspec = 2,
                 names = c("S", "I"), dimens = Gsize, method="euler")
))


times <- seq(1, 10000, by=1/10)
parameters <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu,D=D,D1=0,D2=0)
print(system.time(
  out4 <- ode.1D(y = initial_values, func = SCIR1D,
                 times = times, parms = parameters, nspec = 2,
                 names = c("S", "I"), dimens = Gsize, method="euler")
))


SCIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -b*S + q(I,b,n,k)*D 
    dI <- (beta*S*I)/D - rmu*I
    dR <- rmu*I
    list(c(dS,dI,dR))
  })
}
t <- seq(1, 100000, by=1/10)
init <- c(S = (N-1)/A, I=1/A, R=0)
params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu,D=D)
fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params,method="rk4"))


plot(out1[,(1)+Gsize%/%2+100],out1[,(Gsize+1)+Gsize%/%2+100],xlab="S",ylab="I",type="l",col="blue")
lines(out2[,(1)+Gsize%/%2+10],out2[,(Gsize+1)+Gsize%/%2+10],xlab="S",ylab="I",type="l",col="red")
lines(out3[,(1)+Gsize%/%2+300],out3[,(Gsize+1)+Gsize%/%2+300],xlab="S",ylab="I",type="l",col="green")
lines(out4[,(1)+Gsize%/%2],out4[,(Gsize+1)+Gsize%/%2],xlab="S",ylab="I",type="l",col="purple")

plot(out1[,(1)+Gsize%/%2],out1[,(Gsize+1)+Gsize%/%2],ylim=c(0,100),xlab="S",ylab="I",type="l",col="blue")
lines(out2[,(1)+Gsize%/%2],out2[,(Gsize+1)+Gsize%/%2],ylim=c(0,100),xlab="S",ylab="I",type="l",col="red")
lines(out3[,(1)+Gsize%/%2],out3[,(Gsize+1)+Gsize%/%2],ylim=c(0,100),xlab="S",ylab="I",type="l",col="green")
lines(out4[,(1)+Gsize%/%2],out4[,(Gsize+1)+Gsize%/%2],ylim=c(0,100),xlab="S",ylab="I",type="l",col="purple")

lines(fit$S,fit$I,xlab="S",ylab="I",type="l",col="purple")





dt1<-data.frame(Diffusion=as.factor(c(rep("yes",length(out1[,(Gsize+1)+Gsize%/%2+100])),rep("no",length(out4[,(Gsize+1)+Gsize%/%2])))),I=c(out1[,(Gsize+1)+Gsize%/%2+100],out4[,(Gsize+1)+Gsize%/%2]),S=c(out1[,(1)+Gsize%/%2+100],out4[,(1)+Gsize%/%2]))

ggplot(dt1, aes(x=I, y=S, color=Diffusion)) +
  geom_path() + 
  theme_minimal() +
  xlab(expression(paste("Infected individuals/", km^2))) + 
  ylab(expression(paste("Susceptible individuals/", km^2))) +
  ggtitle("Simulations with and without diffusion") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 20))


ggplot(dt1, aes(x=I, y=S, color=Diffusion)) +
  geom_path() + 
  theme_minimal() +
  xlim(0,1) +
  ylim(0,7) +
  xlab(expression(paste("Infected individuals/", km^2))) + 
  ylab(expression(paste("Susceptible individuals/", km^2))) +
  ggtitle("Simulations with and without diffusion") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 20))
