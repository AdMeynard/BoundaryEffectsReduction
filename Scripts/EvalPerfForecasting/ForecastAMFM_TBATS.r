library('forecast')

LoadConvertTS <- function(TabName) {
  x = read.table(TabName)
  y = ts(x)
  return(y)
}			

#setwd(getSrcDirectory(LoadConvertTS)[1])

### TBATS Extension

TBATSextension <- function(sig,ss,L,xext0) {
  matExt = matrix(0,nrow=length(ss),ncol=L)
  compVal = rep(0,length(ss))
  k = 1
  for (s in ss) {
    fit <- tbats(sig,use.box.cox=FALSE,use.trend=FALSE,use.damped.trend=FALSE,seasonal.periods=s)
    TBATSest = forecast( fit, h=L )
    xext = as.numeric( TBATSest[["mean"]] )
    matExt[k,] = xext
    compVal[k] = sum((xext - xext0)^2)
    k = k+1
  }
  kopt = which.min(compVal)
  xextOpt = matExt[kopt,]
  
  xxTBATS = c( x, xextOpt)
  return(xxTBATS)
}

# Signal
xx0 = LoadConvertTS("../../Signals/AMFM4R")

N = 10000
fs = N-1
t = seq(from = 0, to = 1, length.out = N)

### forecasting parameters
HOP = 1
extSEC = 0.1 # the extension is of extSEC second
L = round( extSEC*fs )
xx0 = xx0[(L+1):(N+2*L)] # remove the left boundary/extension
tt = seq(from = 0, to = 1+L/fs, by = 1/fs)

x0 = xx0[1:N] # restriction to the measurement interval
xextR = xx0[(N+1):(N+L)]

sigman = 1e-2

# Forecasting
nbXP = 15
ss = 1:250

TBATStime = rep(0,nbXP)
MeanTBATS = matrix(0,nrow=nbXP,ncol=L)
VarTBATS = matrix(0,nrow=nbXP,ncol=L)

for (ind in 1:nbXP){
  noise = sigman*rnorm(N+L)
  x = x0 + noise[1:N] # signal to be extended
  
  start_time = Sys.time()
  xxTBATS = TBATSextension(x,ss,L,xextR)
  end_time = Sys.time()
  
  TBATStime[ind] = 2*(end_time - start_time)
  MeanTBATS[ind,] = xxTBATS[(N+1):(N+L)] - xx0[(N+1):(N+L)]
  VarTBATS[ind,] = ( xxTBATS[(N+1):(N+L)] - xx0[(N+1):(N+L)] )^2
}

BiasXP_TBATS = colMeans(MeanTBATS)
VarianceXP_TBATS = colMeans(VarTBATS)
CPUtimeXP_TBATS = mean(TBATStime)

MSE_TBATS = BiasXP_TBATS^2 + VarianceXP_TBATS

length(CPUtimeXP_TBATS) = length(BiasXP_TBATS)
ResTBATS=cbind(BiasXP_TBATS,VarianceXP_TBATS,CPUtimeXP_TBATS)
write.csv(ResTBATS,'../../Results/PerfAHM_TBATS.csv')
