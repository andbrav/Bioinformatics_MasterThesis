
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
b=0.03
k= 150000
rmu=0.021
n=10

t <- seq(1, 100000, by=1/10)
params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu)
init <- c(S = N-1, C=0, I=1, R=0)
fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params,method="rk4"))

print(ggplot(fit, aes(x=I, y=S)) +
geom_hline(aes(linetype = "I Nullcline", yintercept=rmu*N/beta),color="blue") +
geom_vline(xintercept = 0, color="blue",lty=2) +
geom_path(aes(color=C),alpha=1) + 
scale_colour_gradient(name="Confined",low = "green", high = "red") +
geom_vline(aes(linetype = "K",xintercept = k), show.legend = F) +
geom_vline(aes(linetype = "Prediction",xintercept = k*((beta-rmu)/rmu)^(1/n)), show.legend = F, color="purple") +
xlab("Number of infected people") + 
ylab("Number of susceptible people") +
xlim(0,400000) +
ylim(0,5000000) +
theme_minimal() +
ggtitle("Numerical Integration with Confinement Policy") +
labs(subtitle="n=10") +
theme(plot.title = element_text(hjust = 0.5)) +
theme(plot.subtitle = element_text(hjust = 0.5)) +
scale_linetype_manual(name = "Legend", values = c(2,1,1), 
   guide = guide_legend(order=2,override.aes = list(color = c("blue","black" ,"purple")))))+
theme(text = element_text(size = 20)) 



t <- seq(1, 100000, by=1/10)
#b=0 so there is no confinement
params <- c( N=N,beta=beta,b=0,n=n,k=k,rmu=rmu)
init <- c(S = N-1, C=0, I=1, R=0)
fit2 <- data.frame(ode(y = init, times = t, func = SCIR, parms = params,method="rk4"))
print(ggplot(fit2, aes(x=I, y=S)) +
        geom_path(aes(color=C),alpha=1) + 
        scale_colour_gradient(name="Confined", limits=c(0,16000000) ,low = "green", high = "red", breaks=c(0,4000000,8000000,12000000,16000000)) +
        geom_hline(aes(linetype = "I Nullcline", yintercept=rmu*N/beta),color="blue") +
        geom_vline(aes(linetype = "K",xintercept = k), show.legend = F) +
        geom_vline(xintercept = 0, color="blue",lty=2) +
        xlab("Number of infected people") + 
        ylab("Number of susceptible people") +
        xlim(0,400000) +
        ylim(0,5000000) +
        theme_minimal() +
        ggtitle("Numerical Integration without Confinement Policy") +
        #labs(subtitle="n=30") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.subtitle = element_text(hjust = 0.5)) +
        scale_linetype_manual(name = "Legend", values = c(2,1), 
                              guide = guide_legend(order=2,override.aes = list(color = c("blue","black")))))+
  theme(text = element_text(size = 20)) 

# Time plot

nlist=seq(1,40,1)
blist=seq(0,0.12,0.01)
#The second condition is to avoid getting the start of the epidemic
times=c(fit[which(fit2$I<100000 & fit2$S< 20000000),][1,1])
df <- data.frame(b=double(),n=double(),time=double())

for(n in nlist){
  for(b in blist){
    
       t <- seq(1, 100000, by=1/10)
       params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu)
       init <- c(S = N-1, C=0, I=1, R=0)
       fit3 <- data.frame(ode(y = init, times = t, func = SCIR, parms = params,method="rk4"))
       row <- data.frame(b=b,
                         n=n,
                         time=fit[which(fit3$I<100000 & fit3$S< 20000000),][1,1])
       print(row)
       df <- rbind(df,row)
  }
}

saveRDS(df,"dff.rds")
df<-readRDS("dff.rds")
noconft <- fit[which(fit2$I<100000 & fit2$S< 20000000),][1,1]
data_times<-data.frame(x=c(0,nlist),y=times)
dt0<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0,]$time))
colnames(dt0)<-c("n","Days")
dt1<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.01,]$time))
colnames(dt1)<-c("n","Days")
dt2<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.02,]$time))
colnames(dt2)<-c("n","Days")
dt3<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.03,]$time))
colnames(dt3)<-c("n","Days")
dt4<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.04,]$time))
colnames(dt4)<-c("n","Days")
dt5<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.05,]$time))
colnames(dt5)<-c("n","Days")
dt6<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.06,]$time))
colnames(dt6)<-c("n","Days")
dt7<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.07,]$time))
colnames(dt7)<-c("n","Days")
dt8<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.08,]$time))
colnames(dt8)<-c("n","Days")
dt9<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.09,]$time))
colnames(dt9)<-c("n","Days")
dt10<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.10,]$time))
colnames(dt10)<-c("n","Days")
dt11<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.11,]$time))
colnames(dt11)<-c("n","Days")
dt12<-data.frame(x=c(0,nlist),y=c(noconft,df[df$b==0.12,]$time))
colnames(dt12)<-c("n","Days")

data <- cbind(rbind(dt3,dt6,dt9,dt12),as.factor(c(rep(0.03,41),rep(0.06,41),rep(0.09,41),rep(0.12,41))))
colnames(data)<-c("n","Days","b")
ggplot(data=data, aes(x=n, y=Days,color=b))+
        geom_point()+
        theme_minimal() +
        ggtitle("Time until Infected cases < 100000") +
        ylab("Time") +
        theme(plot.subtitle = element_text(hjust = 0.5)) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        scale_y_continuous(limits=c(0,12500),breaks = c(0,1092,2184,3276,4368,5460,6552,7644,8736,9828,10920,12012,13104,14196),
               labels=c("0"="0 yr","1092" = "3 yr","2184" = "6 yr","3276" = "9 yr","4368" = "12 yr","5460" = "15 yr","6552"="18 yr","7644"="21 yr","8736"="24 yr","9828"="27 yr","10920"="30 yr","12012"="33 yr","13104"="36 yr","14196"="39 yr"))+
        scale_x_continuous(breaks = c(0, 10, 20,30,40),labels=c("0" = "No Confinement","10"="10","20"="20","30"="30","40"="40"))+
  xlab("n (harshness of change)")+
  scale_color_discrete(name="b (speed)")+
  theme(text = element_text(size = 20)) 

