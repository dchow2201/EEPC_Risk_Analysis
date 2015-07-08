EEPC <- read.csv("../Data/EEPC Risk Analysis.csv", stringsAsFactors=FALSE)
Total <- read.csv("../Data/EEPC Total.csv", stringsAsFactors=FALSE)
set.seed(123)
summary(EEPC$Amount)
library(moments)
library(dplyr)
library(GLDEX)
library(MASS)
skewness(EEPC$Amount)
kurtosis(EEPC$Amount)

#Exponential distribution
N=nrow(EEPC)
exp_rate=fitdistr(EEPC$Amount,densfun = "exponential")
exp_rate$estimate
exp_ml=1/mean(EEPC$Amount)
hist(EEPC$Amount,xlim = c(0,1e6),breaks = 5e2,freq = F,
     main="Histogram with Exponential Distribution",
     col="green",xaxt="n",yaxt="n",
     xlab="EEPC Law Suit Cost (in thousands)")
xrange=seq.int(0,1e6,by = 100)
yexp=dexp(xrange,exp_rate$estimate)
lines(x = xrange,y = yexp,lwd=3,col="blue")
xrange=seq.int(0,1e6,by = 2e5)
xlabel=paste("$",xrange/1000)
axis(1,at=xrange,labels = xlabel)
yrange=seq.int(0,1.2e-5,by = 2e-6)
ylabel=paste(yrange*100,"%",sep="")
axis(2,at=yrange,labels = ylabel)

rand_exp=rexp(1e5,rate = exp_rate$estimate)
exprange=c(0,1e6)
qqplot(EEPC$Amount,rand_exp,col="darkgrey",pch=19,xlim=exprange,ylim=exprange,
       ylab="Exponential Distribution",xlab="EEPC Law Suit Costs",main="QQ-plot")
abline(0,1,col="red",lty=3,lwd=2)

ks.gof(EEPC$Amount,"pexp",rate = exp_rate$estimate)


#LogNormal
EEPC$LogAmount=log(EEPC$Amount)
lnorm_rate=fitdistr(EEPC$Amount,densfun = "log-normal")
lnorm_rate$estimate
hist(EEPC$Amount,xlim = c(0,1e6),breaks = 5e2,freq = F,
     main="Histogram with Log-Normal Distribution",
     col="yellow",xaxt="n",yaxt="n",
     xlab="EEPC Law Suit Cost (in thousands)")
xrange=seq.int(0,1e6,by = 100)
ylog=dlnorm(xrange,meanlog = lnorm_rate$estimate["meanlog"],
            sdlog = lnorm_rate$estimate["sdlog"])
lines(x = xrange,y = ylog,lwd=3,col="purple")
xrange=seq.int(0,1e6,by = 2e5)
xlabel=paste("$",xrange/1000)
axis(1,at=xrange,labels = xlabel)
yrange=seq.int(0,1.4e-5,by = 2e-6)
ylabel=paste(yrange*100,"%",sep="")
axis(2,at=yrange,labels = ylabel)


ks.gof(EEPC$Amount,"plnorm",meanlog = lnorm_rate$estimate["meanlog"],
       sdlog = lnorm_rate$estimate["sdlog"])


norm_rate=fitdistr(EEPC$LogAmount,densfun = "normal")
norm_rate$estimate
EEPC %>% summarise(mean=mean(LogAmount),sd=sd(LogAmount))
hist(EEPC$LogAmount,freq = F,xlim = c(6,14),breaks = 20,
     main="Histogram with Normal Distribution",
     col="pink",yaxt="n",
     xlab="EEPC Law Suit Cost (Log Scale)")
xrange=seq.int(0,18,by = .1)
ynorm=dnorm(xrange,mean = norm_rate$estimate["mean"],
            sd = norm_rate$estimate["sd"])
lines(x = xrange,y = ynorm,lwd=3,col="darkgreen")

yrange=seq.int(0,0.4,by = 0.05)
ylabel=paste(yrange*100,"%",sep="")
axis(2,at=yrange,labels = ylabel,las=2,cex=.5)

rand_lnorm=rlnorm(1e5,meanlog = lnorm_rate$estimate["meanlog"],
                 sdlog = lnorm_rate$estimate["sdlog"])
exprange=c(0,1e6)
qqplot(EEPC$Amount,rand_lnorm,
       col="darkblue",pch=19,xlim=exprange,ylim=exprange)
abline(0,1,col="red",lty=3,lwd=3)

#Loss Frequence Model
year=EEPC %>% group_by(year.settled) %>% summarise(AvgAmount=mean(Amount),N=n())
library(moments)
year %>% summarise(Avg=mean(N),Var=var(N),Median=median(N),
                   skewness(N),kurtosis(N))
#Because sample variance is much larger than sample mean
#use Negative Binomial


nb_param=fitdistr(year$N,"Negative Binomial")
rnbinom(n=8,size = nb_param$estimate["size"],mu = nb_param$estimate["mu"])

ks.gof(year$N,"pnbinom",size = nb_param$estimate["size"],
       mu = nb_param$estimate["mu"])




#Aggregated loss
total_cost=rep(0,5e4)
for(i in 1:5e4){
    num_case=rnbinom(n=1,size = nb_param$estimate["size"],
                     mu = nb_param$estimate["mu"])
    for(j in 1:num_case){
        total_cost[i]=total_cost[i]+
            rlnorm(n=1,meanlog = lnorm_rate$estimate["meanlog"],
                   sdlog = lnorm_rate$estimate["sdlog"])
    }
}
summary(total_cost)
lnorm_parm=fitdistr(total_cost,densfun = "log-normal")
hist(total_cost,xlim = c(0,2e7),ylim=c(0,2.5e-7),breaks = 1e2,freq = F,
     main="Histogram with Log-Normal Distribution",
     col="yellow",xaxt="n",#yaxt="n",
     xlab="EEPC Law Suit Cost (in millions)")
xrange=seq.int(0,2e7,by = 100)
ylog=dlnorm(xrange,meanlog = lnorm_parm$estimate["meanlog"],
            sdlog = lnorm_parm$estimate["sdlog"])
lines(x = xrange,y = ylog,lwd=3,col="purple")
xrange=seq.int(0,2e7,by = 5e6)
xlabel=paste("$",xrange/1e6,sep="")
axis(1,at=xrange,labels = xlabel)


Total[,c(1,4)]
sumary(Total$Total)

#Kolmogorov-Smirnov test
library(GLDEX)
ks.gof(Total$Total,total_cost)
ks.gof(Total$Total,"plnorm",meanlog = lnorm_parm$estimate["meanlog"],
       sdlog = lnorm_parm$estimate["sdlog"])



library(ggplot2)

ggplot(EEPC,aes(x=LogAmount))+
    geom_histogram(aes(y=..density..),fill="blue",alpha=0.5,binwidth=1)+
    geom_density(color="yellow")+xlab(expression(bold("Law Suit Cost (in log scale)")))+
    ylab(expression(bold("Density")))


#Gamma Distribution
gamma_param=EEPC %>% summarise(scale=var(Amount)/mean(Amount),
                               shape=(mean(Amount)/sd(Amount))^2)

# gamma_rate=fitdistr(EEPC$Amount,dgamma, list(scale = 1e7,shape = 1e-2),
#                     lower=1e-5)
# 
# gamma_rate$estimate

hist(EEPC$Amount,xlim = c(0,1e6),breaks = 5e2,freq = F,
     main="Histogram with Gamma Distribution",
     col="palegreen",xaxt="n",yaxt="n",
     xlab="EEPC Law Suit Cost (in thousands)")
xrange=seq.int(0,1e6,by = 100)
ylog=dgamma(xrange,shape = gamma_param$shape,scale = gamma_param$scale)
lines(x = xrange,y = ylog,lwd=3,col="violetred")
xrange=seq.int(0,1e6,by = 2e5)
xlabel=paste("$",xrange/1000)
axis(1,at=xrange,labels = xlabel)
yrange=seq.int(0,1.4e-5,by = 2e-6)
ylabel=paste(yrange*100,"%")
axis(2,at=yrange,labels = ylabel)

#Weibull distribution
wb_param=fitdistr(EEPC$Amount,densfun = "weibull")

tila_X=function(n){
    sum(log(log(1/(1-(1:n)/(n+1)))))/n
}





n=nrow(EEPC)
EEPC=EEPC %>% mutate(P=rank(Amount,ties.method = "min"))


wb_beta=function(){
    
}








hist(EEPC$Amount,xlim = c(0,1e6),breaks = 5e2,freq = F,
     main="Histogram with Weibull Distribution",
     col="green",xaxt="n",yaxt="n",
     xlab="EEPC Law Suit Cost (in thousands)")
xrange=seq.int(0,1e6,by = 100)
yweibull=dweibull(xrange,shape = wb_param$estimate["shape"],
                  scale = wb_param$estimate["scale"])
lines(x = xrange,y = yweibull,lwd=3,col="blue")
xrange=seq.int(0,1e6,by = 2e5)
xlabel=paste("$",xrange/1000)
axis(1,at=xrange,labels = xlabel)
yrange=seq.int(0,1.2e-5,by = 2e-6)
ylabel=paste(yrange*100,"%",sep="")
axis(2,at=yrange,labels = ylabel)

rand_weibull=rweibull(1e4,shape = wb_param$estimate["shape"],
                  scale = wb_param$estimate["scale"])
exprange=c(0,1e6)
qqplot(EEPC$Amount,rand_weibull,
       col="turquoise",pch=19,xlim=exprange,ylim=exprange)







EEPC$zAmount=scale(EEPC$Amount)


estBetaParams <- function(x) {
  mu=mean(x)
  var=var(x)
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}



logNormParams=fitdistr(EEPC$Amount,densfun = "log-normal")



qqplot(EEPC$Amount,rlnorm(N,meanlog = ))


betaParams=estBetaParams(EEPC$LogAmount)
qqplot(EEPC$LogAmount,rbeta(N,betaParams$alpha,betaParams$beta))

summary(EEPC$LogAmount)
hist(EEPC$LogAmount)