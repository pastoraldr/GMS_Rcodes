rm(list=ls())
library(optimx)
library(signal)
library(readxl)
library(stringr)
library(numbers)
library(latticeExtra)
source("RungeKutta.R")
source("Simulation_Models")


Data <- read_xls("XuRenHua.xls")
colnames(Data)[c(1,2)]=c("time", "glucose")
date <-  NULL;
hr <- NULL;
minute <- NULL;
for(i in 1:dim(Data)[1]){
  tmp <- gsub("[: -]", "" , Data$time[i], perl=TRUE)
  tmp1 <- as.numeric(substr(tmp,5, 8))
  tmp2 <- substr(tmp,9,10)
  tmp3 <- substr(tmp,11,12)
  date = c(date, tmp1)
  hr = c(hr, tmp2)
  minute = c(minute, tmp3)
}
v <- format(as.POSIXct(Data$time,format='%Y-%m-%d %H:%M'),format='%Y/%m/%d')
v = unique(v)

diffDate = as.Date(v, format="%Y/%m/%d")-as.Date(v[1], format="%Y/%m/%d")

Data= data.frame(glucose=Data$glucose, date=as.numeric(date),  hour = as.numeric(hr), minute = as.numeric(minute))
date = unique(Data$date)
overalltime <- rep(0, dim(Data)[1])
for(i in 1:length(date)){
  pos=which(Data$date==date[i]);
  overalltime[pos]=diffDate[i]*1440+Data$hour[pos]*60+Data$minute[pos]
}
Data$overalltime = overalltime

approCurve <- approx(Data$overalltime, y = Data$glucose, xout=seq(min(Data$overalltime), max(Data$overalltime),
                                                                  by=1 ),
                     method="linear",rule = 1, f = 0, ties = mean)

#############################
#############################
#############################
#############################
#############################

NewData <- approCurve
NewData$y = NewData$y*18
idx = which(NewData$x>=1440 & NewData$x<=1440+60*3)
PartNewData <- data.frame(t = NewData$x[idx]-1440, u = NewData$y[idx])


idx = which(NewData$x==1440)
tmp <-  scan("tmp.txt");
tmphistory <-  c( c((1440-60):1440),  
               NewData$y[c((idx-60):idx)],
               tmp[c((183-60):183)])
tspan1 = c(1440, 1440+60*3)
u1 = tmphistory[c(61*2,61*3)]


lower = rep(0, 12)
upper = rep(0, 12)
VP.di = 0.02;   lower[1] = 0.01; upper[1] = 0.08;
VP.sig1 = 10;   lower[2] = 12;   upper[2] = 72;
VP.a1 = 160;    lower[3] = 150;  upper[3] = 210
VP.r1 = 6;      lower[4] = 3;    upper[4] = 8; 
VP.sig2 = 0.64; lower[5] = 0.1;  upper[5] = 0.8;
VP.a2 = 9;      lower[6] = 5;    upper[6] = 50;
VP.sig4 = 0.0099; lower[7] = 1e-6;   upper[7] = 2e-2;
VP.U0 = 0;        lower[8] = -1;     upper[8] = 0.4;
VP.a4 = 50.6;     lower[9] = 0;      upper[9] = 60;
VP.sig5 = 2.25;   lower[10] = 0;     upper[10] = 500;
VP.a5 = 26;       lower[11] = 5;     upper[11] =50;
VP.r5 = 8;        lower[12] = 2;     upper[12] = 9;

lower = c(lower)
upper = c(upper)
v= c(0.25, 40 , 60, 1, 1, 1, 0.004, 1, 20, 1, 20, 5)

tv1 = 1e+5;
result=NULL;
for(k in 15:60){
  tsf <- function(vx){
    vpx = c(2,k,2,vx, 2e+6, tmphistory);
    tmp_sol = rk_fixed(f,u1, h, p=vpx,tspan1, k=1);
    return(mean((tmp_sol$u[,1]- PartNewData$u)^2));
  }
  
  llxx = nlminb(v, tsf, lower= lower, upper=upper, control = list(iter.max = 40))
  result=c(result,llxx$objective)
}

which.min(result)
tsf <- function(vx){
  vpx = c(2,14+8,2,vx, 2e+6, tmphistory);
  tmp_sol = rk_fixed(f,u1, h, p=vpx,tspan1, k=1);
  return(mean((tmp_sol$u[,1]- PartNewData$u)^2));
}

llxx = nlminb(v, tsf, lower= lower, upper=upper, control = list(iter.max = 40))
vpx = c(2,14+46,2,llxx$par, 2e+6, tmphistory);
tmp_sol = rk_fixed(f,u1, h, p=vpx,c(1440,1440+60*5), k=1);

par(mfrow=c(1,1))
idx = which(NewData$x>=1440 & NewData$x<=1440+60*5)
NewData2 <- data.frame(t = NewData$x[idx]-1440, u = NewData$y[idx])
plot(NewData2$t, NewData2$u, type="l")
lines(tmp_sol$t-1440, tmp_sol$u[,1],col=2)

