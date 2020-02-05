function xx = SigExtension(x,fs,HOP,extN,extP,extSEC,method)

xext = forecasting(x,fs,HOP,extN,extP,extSEC,method);
xexti = forecasting(x,fs,HOP,extN,extP,extSEC,method,'backward');

xx = [xexti; x; xext];