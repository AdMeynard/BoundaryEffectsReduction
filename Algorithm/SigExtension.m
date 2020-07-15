function xx = SigExtension(x,fs,HOP,extK,extM,extSEC,method)

L = round(fs*extSEC);

xext = forecasting(x,L,HOP,extK,extM,method);
xexti = forecasting(x,L,HOP,extK,extM,method,'backward');

xx = [xexti; x; xext];