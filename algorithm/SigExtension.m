function xx = SigExtension(x,fs,HOP,extN,extP,extSEC)

xext = forecasting(x,fs,HOP,extN,extP,extSEC);
xexti = forecasting(x,fs,HOP,extN,extP,extSEC,'backward');

xx = [xexti; x; xext];