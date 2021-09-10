
require(deSolve)
require(ggplot2)

p <- function(X,b1,n1,k1) {b1*X**n1/(k1**n1+X**n1)}
q <- function(X,b2,n2,k2) {b2*k2**n2/(k2**n2+ X**n2)}



# Model 1 p(I) and q(I)



SCIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -(beta*S*I)/N + q(I,b,n,k)*C - p(I,b,n,k)*S
    dC <- p(I,b,n,k)*S -  q(I,b,n,k)*C
    dI <- (beta*S*I)/N - rmu*I
    dR <- rmu*I
    list(c(dS,dC,dI,dR))
  })
}

N=45000000
beta=0.45
k= 150000
rmu=0.021
n=10
# Infected plot

nlist=seq(1,40,1)
blist=seq(0.01,0.12,0.01)
df1 <- data.frame(b=double(),n=double(),recovered=double(),max_infected=double())
for(n in nlist){
  for(b in blist){
print(n)
print(b)
t <- seq(1, 100000, by=1/10)
params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu)
init <- c(S = N-1, C=0, I=1, R=0)
fit3 <- data.frame(ode(y = init, times = t, func = SCIR, parms = params,method="rk4"))
row <- data.frame(b=b,
                  n=n,
                  recovered=fit3$R[length(fit3$R)],
                  max_infected = max(fit3$I))
print(row)
df1 <- rbind(df1,row)
}
}

saveRDS(df1,"df1.rds")
df<-readRDS("df1.rds")
t <- seq(1, 100000, by=1/10)
init <- c(S = N-1, C=0, I=1, R=0)
params <- c( N=N,beta=beta,b=0,n=n,k=k,rmu=rmu)
fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params,method="rk4"))
data <- cbind(rbind(dt3,dt6,dt9,dt12),as.factor(c(rep(0.03,41),rep(0.06,41),rep(0.09,41),rep(0.12,41))))


elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})

dt1<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.01),]$max_infected))
dt2<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.02),]$max_infected))
dt3<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.03),]$max_infected))
dt4<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.04),]$max_infected))
dt5<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.05),]$max_infected))
dt6<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.06),]$max_infected))
dt7<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.07),]$max_infected))
dt8<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.08),]$max_infected))
dt9<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.09),]$max_infected))
dt10<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.10),]$max_infected))
dt11<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.11),]$max_infected))
dt12<-data.frame(x=seq(0,40,1),y=c(max(fit$I),df[elementwise.all.equal(df$b,0.12),]$max_infected))
data <- cbind(rbind(dt1,dt3,dt6,dt9,dt12),as.factor(c(rep(0.01,41),rep(0.03,41),rep(0.06,41),rep(0.09,41),rep(0.12,41))))
data2 <- cbind(rbind(dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9,dt10,dt11,dt12),as.factor(c(rep(0.01,41),rep(0.02,41),rep(0.03,41),rep(0.04,41),rep(0.05,41),rep(0.06,41),rep(0.07,41),rep(0.08,41),rep(0.09,41),rep(0.10,41),rep(0.11,41),rep(0.12,41))))
colnames(data)<-c("n","MaxInfected","b")
colnames(data2)<-c("n","MaxInfected","b")

ggplot(data=data, aes(x=n, y=MaxInfected,color=b))+
  geom_point()+
  theme_minimal() +
  ggtitle("Infected individuals peak") +
  ylab("Infected individuals") +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits=c(0,40000000))+
  scale_x_continuous(breaks = c(0, 10, 20,30,40),labels=c("0" = "No Confinement","10"="10","20"="20","30"="30","40"="40"))+
  xlab("n (harshness of change)")+
  scale_color_discrete(name="b (speed)")+
  theme(text = element_text(size = 20)) 

