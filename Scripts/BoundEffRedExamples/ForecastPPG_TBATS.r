library('forecast')

#x = read.table("sigR")
x = read.table("sigR2")

y = ts(x)
N = length(y)

fit <- tbats(y,use.box.cox=FALSE,use.trend=FALSE,use.damped.trend=FALSE,seasonal.periods=23)
plot( forecast( fit, h=floor(N/3) ) )
