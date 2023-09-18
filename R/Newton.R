newton=function(f,xInit,maxIt=20,relConvCrit=1.e-10,...){
  p=length(xInit)
  results=matrix(NA,maxIt,p+2) 
  colnames(results)=c("value",paste("x",1:p,sep=""),"Conv")
  
  xCurrent=xInit
  for(t in 1:maxIt){
    evalF=f(xCurrent,der=2,...)
    results[t,"value"]=evalF$value
    results[t,1+(1:p)]=xCurrent
    xNext=xCurrent-solve(evalF$der2,evalF$der1)
    conv=sqrt(crossprod(xNext-xCurrent))/(sqrt(crossprod(xCurrent))+relConvCrit)
    results[t,"Conv"]=conv
    if(conv < relConvCrit) break
    xCurrent=xNext
  }
  return(list(x=xNext,value=f(xNext,...),convergence=(Conv < relConvCrit),
              results=results[1:t,]))
}

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

source("R/WolfeFunction.R")
newtonWolfe=function(f,xInit,maxIt=20,relConvCrit=1.e-10,...){
  p=length(xInit)
  results=matrix(NA,maxIt,p+2) 
  colnames(results)=c("value",paste("x",1:p,sep=""),"Conv")
  
  xCurrent=xInit
  for(t in 1:maxIt){
    evalF=f(xCurrent,der=2,...)
    results[t,"value"]=evalF$value
    results[t,1+(1:p)]=xCurrent
    p=-solve(evalF$der2,evalF$der1)
    xNext=Wolfe(f,p,xCurrent,...)
    Conv=sqrt(crossprod(xNext-xCurrent))/(sqrt(crossprod(xCurrent))+relConvCrit)
    results[t,"Conv"]=Conv
    if(Conv < relConvCrit) break
    xCurrent=xNext
  }
  return(list(x=xNext,value=f(xNext,...),convergence=(Conv < relConvCrit),
              results=results[1:t,]))
}