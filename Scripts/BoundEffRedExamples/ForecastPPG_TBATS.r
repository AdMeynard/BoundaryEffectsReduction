library('forecast')

LoadConvertTS <- function(TabName) {
  x = read.table(TabName)
  y = ts(x)
  return(y)
}			

setwd(getSrcDirectory(LoadConvertTS)[1])

y = LoadConvertTS("../../Signals/PPG4R")
N = length(y)

fit <- tbats(y,use.box.cox=FALSE,use.trend=FALSE,use.damped.trend=FALSE,seasonal.periods=23)
plot( forecast( fit, h=floor(N/3) ) )