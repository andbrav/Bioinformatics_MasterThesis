


require(deSolve)
library(ggplot2)
library(plot3D)

p <- function(X,b1,n1,k1) {b1*X**n1/(k1**n1 + X**n1)}
q <- function(X,b2,n2,k2) { b2*k2**n2/(k2**n2 + X**n2)}

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




N = 45000000
k = 150000
beta = 0.45
rmu = 0.021

blist = seq(0.01,0.12,0.01)
nlist = seq(1,30,1)

fs <- 100 # sampling frequency
t <- seq(1, 2000, by=1/fs)
df1 <- data.frame(b=double(),n=double(),oscillations=double(),fourier = double(),
                 prediction2d = double())
df2 <- data.frame(b=double(),n=double(),oscillations=double(),fourier = double(),
                  prediction2d = double())

i=1
signals<-list() # We will store all the integrations to save computation time
padding<-60000*10 # We use padding to get more resolution from fourier analysis

for(b in blist){ 
  for(n in nlist){
    cond <-  (2*beta*(beta/rmu-1))/(beta/rmu)^2
    oscillations <- 0 # oscillation flag (we are only interested in conditions that give oscillatory behavior)
    if(((b<cond) & n>=1)|(b>cond & n> b*(beta/rmu)^2/(2*beta*(beta/rmu-1)))){
      oscillations <- 1
    }
    print(n)
    print(b)
    print(oscillations)
    params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu)
    init <- c(S = N-1,C=0, I=1,R=0)
    fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params))
    original_signal <- fit$S[(301*fs-(fs-1)):(1000*fs -(fs-1))]
    signals[[i]]<-original_signal
    # We remove the trend with two different methods and will keep the one that gives better results for each case
    signal1 <- c(diff(original_signal,lag=1,differences = 1)-mean(diff(original_signal,lag=1,differences = 1)),rep(0,padding))
    signal2 <- c(diff(original_signal,lag=1,differences = 2)-mean(diff(original_signal,lag=1,differences = 2)),rep(0,padding))
    

    X.k <- fft(signal1)
    freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
    X.k <- X.k[1:(length(X.k)/2)]
    fourierf <- freqrad[which.max(Mod(X.k))]
    prediction2d <- sqrt(b)*sqrt(-(b+4*n*(rmu/beta)*(-beta+rmu)))*(1/2)
    row <- data.frame(b=b,
                      n=n,
                      oscillations=oscillations,
                      fourier = fourierf,
                      prediction2d = prediction2d)
    print(row)
    df1 <- rbind(df1,row)
    
    
    X.k <- fft(signal2)
    freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
    X.k <- X.k[1:(length(X.k)/2)]
    fourierf <- freqrad[which.max(Mod(X.k))]
    prediction2d <- sqrt(b)*sqrt(-(b+4*n*(rmu/beta)*(-beta+rmu)))*(1/2)
    row <- data.frame(b=b,
                      n=n,
                      oscillations=oscillations,
                      fourier = fourierf,
                      prediction2d = prediction2d)
    print(row)
    df2 <- rbind(df2,row)
    
    i=i+1
  }
}

df1$prediction2d[which(!is.na(df1$prediction2d))]=0
df2$prediction2d[which(!is.na(df2$prediction2d))]=0


saveRDS(signals,file="signals.rds")
saveRDS(df1,file="dfdiff2.rds")
saveRDS(df2,file="dfdiff1.rds")


signals<-readRDS("signals.rds")
df1<-readRDS("dfdiff1.rds")
df2<-readRDS("dfdiff2.rds")

dfcombined<-data.frame(b=double(),n=double(),oscillations=double(),fourier = double(),
                       prediction2d = double(),origin=factor())
condition <- abs(df2$prediction2d-df2$fourier) < abs(df1$prediction2d-df1$fourier)
condition[which(is.na(condition))] = TRUE

for(i in 1:360){
  if(condition[i]){
    row <- cbind(df2[i,],origin="df2")
    dfcombined <- rbind(dfcombined,row)
    
  }
  else{
    row <- cbind(df1[i,],origin="df1")
    dfcombined <- rbind(dfcombined,row)
    
  }
  
}

# Some of the results are bad because the signal wasnt cut properly or the trend wasnt removed properly.
# We repeat these cases manualy
# b=0.06 n=1 index=151

b=0.06
n=1
fs <- 100
params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu)
init <- c(S = N-1,C=0, I=1,R=0)
fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params))
original_signal <- fit$S[(55*fs-(fs-1)):(1000*fs -(fs-1))]
padding<-60000*10
signal <- diff(original_signal,lag=1,differences = 1)
signal<-signal-mean(signal)
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[151,]$fourier<-fourierf

# b=0.06 n=2 index=152


padding<-60000
fs <- 100
signal <- diff(signals[[152]][4000:length(signals[[152]])],lag=1,differences = 2)
signal<-signal-mean(signal)
time <- seq(1,length(signal),1)
trend <- lm(signal~1+time)

signal <- signal-predict(trend)
signal[which(is.na(signal))]=0

signal<-c(signal,rep(0,padding))
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[152,]$fourier<-fourierf



# b=0.07 n=1 index=181

padding<-60000*10
fs <- 100
signal <- diff(signals[[181]][12000:length(signals[[181]])],lag=1,differences = 1)
time <- seq(1,length(signal),1)
trend <- lm(signal~1+time)

signal <- signal-predict(trend)
signal[which(is.na(signal))]=0
signal<-signal-mean(signal)
signal<-c(signal,rep(0,padding))
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[181,]$fourier<-fourierf

# b=0.07 n=2 index=182

padding<-60000*10
fs <- 100
signal <- diff(signals[[182]][2000:length(signals[[182]])],lag=1,differences = 2)
signal[which(is.na(signal))]=0
signal<-signal-mean(signal)
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[182,]$fourier<-fourierf

# b=0.09 n=3 index=243

padding<-60000*10
fs <- 100
signal <- diff(signals[[243]][2700:length(signals[[243]])],lag=1,differences = 2)
signal[which(is.na(signal))]=0
signal<-signal-mean(signal)
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[243,]$fourier<-fourierf

# b=0.11 n=2 index=302

b=0.11 
n=2
params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu)
init <- c(S = N-1,C=0, I=1,R=0)
fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params))
original_signal <- fit$S[(50*fs-(fs-1)):(1000*fs -(fs-1))]
padding<-60000*10
fs <- 100
signal <- diff(original_signal,lag=1,differences = 1)
time <- seq(1,length(signal),1)
trend <- lm(signal~1+time)

signal <- signal-predict(trend)
signal[which(is.na(signal))]=0
signal<-signal-mean(signal)
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[302,]$fourier<-fourierf


# b=0.11 n=4 index=304

fs <- 100
signal <- diff(signals[[304]][2000:length(signals[[304]])],lag=1,differences = 1)
signal<-signal-mean(signal)
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[304,]$fourier<-fourierf


# b=0.12 n=2 index=332
dfcombined[332,]

N = 45000000
k = 150000
b = 0.03
beta = 0.45
rmu = 0.021
b=0.12 
n=2
params <- c( N=N,beta=beta,b=b,n=n,k=k,rmu=rmu)
init <- c(S = N-1,C=0, I=1,R=0)
fit <- data.frame(ode(y = init, times = t, func = SCIR, parms = params))
fs <- 100
original_signal <- fit$S[(50*fs-(fs-1)):(1000*fs -(fs-1))]
signal <- diff(original_signal,lag=1,differences = 1)
signal<-signal-mean(signal)
X.k <- fft(signal)
freqrad = (seq(0,fs,fs/length(X.k))*2*pi)
X.k <- X.k[1:(length(X.k)/2)]
fourierf <- freqrad[which.max(Mod(X.k))]
dfcombined[332,]$fourier<-fourierf

saveRDS(dfcombined,file="dfcombined.rds")

#### plots

dfcombined<-readRDS("dfcombined.rds")

dfosc<-dfcombined[!is.na(dfcombined$prediction2d),]

N = 45000000
k = 150000
beta = 0.45
rmu = 0.021
x = seq(0.01,0.12,0.01)
y = seq(1,30,1)
z = matrix(data=NA, nrow=length(x), ncol=length(y))
for(i in 1:length(x))
{
  for(j in 1:length(y))
  {  
    z[i,j] = log(1/(sqrt(x[i])*sqrt(-(x[i]+4*y[j]*(rmu/beta)*(-beta+rmu)))*(1/2)))
  }
}
#z[which(is.na(z))]=0

scatter3D(dfosc$b,dfosc$n,log(1/dfosc$fourier),col="black",pch = 18,xlab="b (speed)",ylab="n (starkness)",zlab="Period", ticktype = "detailed",bty="b2",nticks=40,zlim=c(0,7),xlim=c(0,0.13),ylim=c(0,31),surf=list(x = x, y = y, z = z, fit = log(1/dfosc$prediction2d),alpha=0.2),theta=180)



ggplot(dfosc, aes(x=as.factor(n), y=abs(1/fourier-1/prediction2d)/(1/fourier))) + 
  geom_boxplot(aes(fill="red")) +
  geom_text(aes(label=b),  size=5) +
  ylab("Period Error/ \"Reality\"")+
  xlab("n (starkness)") + 
  ggtitle("Relative error size in function of n")+
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),limits=c(0, 0.75))+
  theme_minimal()+
  theme(legend.position = "None",plot.title = element_text(hjust = 0.5),text = element_text(size = 20))


ggplot(dfosc, aes(x=as.factor(b), y=abs(1/fourier-1/prediction2d)/(1/fourier))) + 
  geom_boxplot(aes(fill="red")) +
  geom_text(aes(label=n),  size=6)+
  ylab("Period Error/ \"Reality\"")+
  xlab("b (speed)") + 
  ggtitle("Relative error size in function of b")+
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),limits=c(0, 0.75))+
  theme_minimal()+
  theme(legend.position = "None",plot.title = element_text(hjust = 0.5),text = element_text(size = 20))
