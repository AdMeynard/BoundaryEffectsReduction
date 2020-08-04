library('forecast')

LoadConvertTS <- function(TabName) {
  x = read.table(TabName)
  y = ts(x)
  return(y)
}			

setwd(getSrcDirectory(LoadConvertTS)[1])

x = LoadConvertTS("../../Signals/ECG4R")
xext0 = LoadConvertTS("../../Signals/ECGext4R")

fs = 200
extSEC = 6 # the extension is of extSEC second
L = round( extSEC*fs )

N = length(x)
t = seq(from = 0, to = (N-1)*fs, by = fs)
tt = seq(from = 0, to = (N+L-1)*fs, by = fs)

xn = x+rnorm(N,sd=0.3) # regularization via awgn

ss = seq(from=160, to=195, by=1)
matExt = matrix(0,nrow=length(ss),ncol=L)
compVal = rep(0,length(ss))
k = 1
for (s in ss) {
  fit <- tbats(xn,use.box.cox=FALSE,use.trend=FALSE,use.damped.trend=FALSE,seasonal.periods=s)
  TBATSest = forecast( fit, h=L )
  xext = as.numeric( TBATSest[["mean"]] )
  matExt[k,] = xext
  compVal[k] = sum((xext - xext0)^2)
  k = k+1
}

kopt = which.min(compVal)
xextOpt = matExt[kopt,]

xx0 = c( x, xext0)
xxTBATS = c( x, xextOpt)
plot(tt,xx0,type="l", col="blue", lty=2)
lines(tt,xxTBATS,type="l", col="red", lty=1)
lines(t,x,type="l", col="blue", lty=1)

write.csv(xextOpt,'../../Results/ECG_TBATS.csv')