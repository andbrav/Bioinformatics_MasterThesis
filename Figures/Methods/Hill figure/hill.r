library(ggplot2)
hill1 <- function(x,b,k,n){b*(x^n)/(k^n+x^n)}


ggplot(data = data.frame(x = 0), mapping = aes(x = x))+
  stat_function(fun = hill1, args = list(b = 2, k = 5, n=2),aes(colour = "2")) +
  stat_function(fun = hill1, args = list(b = 2, k = 5, n=5),aes(colour = "5")) +
  stat_function(fun = hill1, args = list(b = 2, k = 5, n=10),aes(colour = "10"))+
  stat_function(fun = hill1, args = list(b = 2, k = 5, n=15),aes(colour = "15")) +
  scale_colour_manual("n", breaks = c("2","5","10","15"), values = c("green", "blue", "red", "black"))+
  geom_hline(yintercept=2, linetype="dashed", color = "orange")+
  geom_vline(xintercept = 5, linetype="dashed", color = "orange")+
  theme_minimal() +
  scale_x_continuous(breaks=c(0,5,10),limits = c(0, 10),labels=c("","k",""))+
  scale_y_continuous(breaks=c(0,1,2),labels=c("","","b"))+
  labs(x ="x",y="f(x)")+
  theme(text = element_text(size = 20)) 
  


hill2 <- function(x,b,k,n){b*(k^n)/(k^n+x^n)}


ggplot(data = data.frame(x = 0), mapping = aes(x = x))+
  stat_function(fun = hill2, args = list(b = 2, k = 5, n=2),aes(colour = "2")) +
  stat_function(fun = hill2, args = list(b = 2, k = 5, n=5),aes(colour = "5")) +
  stat_function(fun = hill2, args = list(b = 2, k = 5, n=10),aes(colour = "10"))+
  stat_function(fun = hill2, args = list(b = 2, k = 5, n=15),aes(colour = "15")) +
  scale_colour_manual("n", breaks = c("2","5","10","15"), values = c("green", "blue", "red", "black"))+
  geom_hline(yintercept=2, linetype="dashed", color = "orange")+
  geom_vline(xintercept = 5, linetype="dashed", color = "orange")+
  theme_minimal() +
  scale_x_continuous(breaks=c(0,5,10),limits = c(0, 10),labels=c("","k",""))+
  scale_y_continuous(breaks=c(0,1,2),labels=c("","","b"))+
  labs(x ="x",y="f(x)")+
  theme(text = element_text(size = 20)) 