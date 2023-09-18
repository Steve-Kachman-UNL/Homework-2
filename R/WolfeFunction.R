Wolfe<-function(g,p,x0,alphaInit=1,alphaMax=2,c1=1.e-4,c2=.9,...){
  g0=g(x0,1,...)
  phi0=g0$value
  pt=t(p)
  phi0.p=pt%*%g0$der1
  if(phi0.p<0){ #descent direction turn around
    p=-p
    pt=-pt
    phi0.p=-phi0.p
  }
  alpha=alphaInit
  x=x0+alpha*p
  if(phi0.p <= .Machine$double.eps){
    return(x)
  }
  alphaLow=Inf
  alphaHigh=0
  phiHigh=phi0
  while(alpha<=alphaMax ){
    x=x0+alpha*p
    gx=g(x,1,...)
    phi=gx$value
    if(phi<phi0+c1*alpha*phi0.p || phi<phiHigh){  # We jumped to the decreasing part
      alphaLow=alpha
      phiLow=phi
      alpha=(alphaLow+alphaHigh)/2
    }else{  #Cond 1 satisfied 
      phi.p=pt%*%gx$der1
      if(abs(phi.p) < c2*phi0.p || abs(alphaHigh-alphaLow)<1.e-2) return(x)
      #Cond 2 not satisfied
      alphaDiff=alphaHigh-alphaLow
      if(phi.p*(alphaHigh-alphaLow)>0){ # Are we on the down hill side?
        alphaLow=alphaHigh
        phiLow=phiHigh
        alphaHigh=alpha
        phiHigh=phi
        alpha=(alphaLow+alphaHigh)/2
      }
      else{
        alphaHigh=alpha
        phiHigh=phi
        if(alphaLow==Inf){
          alpha=alpha*2
        }
        else{
          alpha=(alphaHigh+alphaLow)/2
        }
      }
    }
  }
  return(x)
}

