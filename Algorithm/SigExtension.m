function xx = SigExtension(x,fs,HOP,extK,extM,extSEC,method)

xext = forecasting(x,fs,HOP,extK,extM,extSEC,method);
xexti = forecasting(x,fs,HOP,extK,extM,extSEC,method,'backward');

xx = [xexti; x; xext];