rm(list=ls())
library(optimx)
library(signal)
source("RungeKutta.R")
source("Simulation_Models")
#############################
#############################
#############################
#############################
#############################



compareTwoSol <- function(sol1, sol2){
  G = sol1$u[,1]
  I = sol1$u[,2]
  len = length(G)
  Gtmp = rep(0, len)
  Itmp = rep(0, len)
  for(k in 1:len){
    idx<-which.min(abs(sol1$t[k]-sol2$t))
    Gtmp[k] = sol2$u[idx,1];
    Itmp[k] = sol2$u[idx,2];
  }
  return(c(mean((Gtmp- G)^2), mean((Itmp- I)^2)))
}

VP.tau1= 5;
VP.tau2 = 36;
VP.tau3 = 3;

VP.di = 0.04;

## parameters for f1
VP.sig1 = 75;
VP.a1 = 203;
VP.r1 = 4; 

## parameters for f2
VP.sig2 = 0.24
VP.a2 = 9;

## parameters for f3 and f4
VP.sig4 = 0.009;
VP.U0 = 0.4;     
VP.a4 = 50.6;     

## parameters for f5
VP.sig5 = 1.77; 
VP.a5 = 26;   
VP.r5 = 8;     

Delay = c(VP.tau1,  VP.tau2,  VP.tau3);
MaxDelay = max(Delay);

StartTime = 0;
VPStartTime = 420
VPEndTime = VPStartTime+60*6;


vpinit = c(VP.tau1, VP.tau2, VP.tau3, VP.di,
           VP.sig1, VP.a1, VP.r1,
           VP.sig2, VP.a2,
           VP.sig4, VP.U0, VP.a4,   
           VP.sig5, VP.a5, VP.r5)
initHistory =initData(60)
vp0 = c(vpinit, 2e+6, initHistory$history)
u0 <- initHistory$u
soltrue = rk_fixed(fgen, u0, h, p=vp0, c(0, VPEndTime)) 
Data = soltrue
idx = which(Data$t>=VPStartTime)
Data.sel = Data[idx,]

tmp= SoltoInit(Data,60,VPStartTime)
vp1 = c(vpinit, 2e+6, tmp$history)
u1 <- tmp$u
tspan1 = c(VPStartTime, VPEndTime )
sol = rk_fixed(f,u1, h, p=vp1, tspan1)


###### find the delays
idx= which(Data$t>80)

Diff1 = Data$u[idx,1]-Data$u[idx-10,1]
Diff2 = Data$u[idx,2]-Data$u[idx-10,2]

ll11=ll22=rep(0,80);
ll12=ll21=rep(0,80);
for(k in 1:80){
  ll11[k]=cor(Diff1, Data$u[idx-10*k,1])
  ll12[k]=cor(Diff1, Data$u[idx-10*k,2])
  ll21[k]=cor(Diff2, Data$u[idx-10*k,1])
  ll22[k]=cor(Diff2, Data$u[idx-10*k,2])
}
tsf <- function(v){
  vpx = c(2,which.max(abs(ll12)),20,v, 2e+6, tmp$history);
  tmp_sol = rk_fixed(f,u1, h, p=vpx,tspan1);
  print(mean((tmp_sol$u[,1]- Data.sel$u[,1])^2)); 
  return(mean((tmp_sol$u[,1]- Data.sel$u[,1])^2+ (tmp_sol$u[,2]- Data.sel$u[,2])^2));
}

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
for(k in 1:30){
  llxx = nlminb(v, tsf, lower= lower, upper=upper, control = list(iter.max = 20))
  tv2 = llxx$objective;
  if(abs(tv1-tv2)<1e-3){
    break;
  }
  tv1 = tv2;
  if (llxx$objective/mean(Data.sel$u[,1])<0.02) {
    break;
  } else{
    v = llxx$par
  }
}





vpx = c(2,33,20,llxx$par, 2e+6, tmp$history);
tmp_sol = rk_fixed(f,u1, h, p=vpx,tspan1)



par(mfrow=c(2,1))
plot(Data.sel$t, Data.sel$u[,1], type="l", col=1)
lines(sol$t, sol$u[,1], col=2)
lines(tmp_sol$t, tmp_sol$u[,1], col=4)


plot(Data.sel$t, Data.sel$u[,2], type="l",col=1)
lines(sol$t, sol$u[,2], col=2)
lines(tmp_sol$t, tmp_sol$u[,2], col=4)


compareTwoSol(Data.sel, sol)
compareTwoSol(Data.sel, tmp_sol)

