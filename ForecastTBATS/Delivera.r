library('forecast')

#x = read.table("sigR")
x = read.table("sigR2")

y = ts(x)
N = length(y)

fit <- tbats(y)
plot( forecast( fit, h=floor(N/3) ) )
