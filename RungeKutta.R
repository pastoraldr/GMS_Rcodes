rk <- function(f, u, h, p, t ,delta){
  k1 <- delta*f(u, h, p, t)
  k2 <- delta* f(u+0.5*k1, h, p, t+0.5*delta)
  k3 <- delta* f(u+0.5*k2, h, p, t+0.5*delta)
  k4 <- delta* f(u+k3, h, p, t+delta)
  nu <- u+(k1 + k2 + k2 + k3 + k3 + k4) / 6
  return(c(t+delta,nu))
}

rk_fixed <- function(f, u, h, p, ts, k=10) {
  
  
  tslot <- seq(ts[1], ts[2], by=1/k)
  sol <- data.frame(t = tslot)
  n <- length(tslot)
  uslot <- matrix(0, nrow=n, ncol=2);
  uslot[1,] = u
  delta <- tslot[2]-tslot[1]
  hnew <- function(subp,subt){
    if(subt>=tslot[1]){
      idx = which.min(abs(tslot-subt))
      return(uslot[idx,])
    } else{
      return(h(p,subt))
    }
  }
  for(i in 2:n) {
    tt <- tslot[i-1]
    uu <- uslot[i-1,]
    uslot[i,] = rk(f, uu, hnew, p, tt ,delta)[-1]
  }
  sol$u <- uslot;
  return(sol)
}

rk_adaptive <- function(f, u, h, p, ts, abserr) {
  tsol = c(ts[1])
  usol = u
  p0 = 4
  
  hnew <- function(subp,subt){
    if(subt>=tsol[1]){
      idx = which.min(abs(tsol-subt))
      return(usol[idx,])
    } else{
      return(h(p,subt))
    }
  }
  
  p21 = 2^p0-1
  tend = ts[2]
  tt = tsol[1]
  uu = usol
  while(tt<tend){
    delta = 0.1
    tmp1 = rk(f, uu, hnew, p, tt ,delta)
    tmp2 = rk(f, uu, hnew, p, tt ,delta/2)
    if ( sum(abs(tmp1[-1]-tmp2[-1]))/p21<abserr){
      while(sum(abs(tmp1[-1]-tmp2[-1]))/p21<abserr){
        tmp2 = tmp1
        delta = 2*delta
        tmp1 = rk(f, uu, hnew, p, tt ,delta)
      }
    } else{
      while(sum(abs(tmp1[-1]-tmp2[-1]))/p21>abserr){
        tmp1 = tmp2
        delta = delta/2
        tmp2 = rk(f, uu, hnew, p, tt ,delta)
      }
    }
    tmp = rk(f, uu, hnew, p, tt ,delta)
    tt = tmp[1]
    uu = tmp[-1]
    tsol = c(tsol, tt)
    usol = rbind(usol, uu)
  }
  return(list(t=tsol, u=usol))
}
