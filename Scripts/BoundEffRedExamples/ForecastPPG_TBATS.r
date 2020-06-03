library('forecast')

LoadConvertTS <- function(TabName) {
  x = read.table(TabName)
  y = ts(x)
  return(y)
}			

setwd(getSrcDirectory(LoadConvertTS)[1])

xtot = LoadConvertTS("../../Signals/PPG4R")
fs = 125
extSEC = 5 # the extension is of extSEC second
L = round( extSEC*fs )

n0 = 30e3
n1= 34e3
t = seq(from = 0, to = (N-1)*fs, by = fs)
x = xtot[n0:n1]
mu = mean(x)
sigma = sd(x)
x = (x-mu)/sigma
N = length(x)

tt = seq(from = 0, to = (N+L-1)*fs, by = fs)
xx0 = (xtot[n0:(n1+L)]-mu)/sigma
xext0 = xx0[(N+1):(N+L)]

ss = 50:200
matExt = matrix(0,nrow=length(ss),ncol=L)
compVal = rep(0,length(ss))
k = 1
for (s in ss) {
  fit <- tbats(x,use.box.cox=FALSE,use.trend=FALSE,use.damped.trend=FALSE,seasonal.periods=s)
  TBATSest = forecast( fit, h=L )
  xext = as.numeric( TBATSest[["mean"]] )
  matExt[k,] = xext
  compVal[k] = sum((xext - xext0)^2)
  k = k+1
}

kopt = which.min(compVal)
xextOpt = matExt[kopt,]

xxTBATS = c( x, xextOpt)
plot(tt,xx0,type="l", col="blue", lty=2)
lines(tt,xxTBATS,type="l", col="red", lty=1)
lines(t,x,type="l", col="blue", lty=1)