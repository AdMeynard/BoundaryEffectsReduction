function xx = SigExtension(x,fs,HOP,extK,extM,extSEC,method)
% SIGEXTENSION Extend signal on both sides by forecasting
% Usage:	xext = forecasting(x,L,HOP,extK,extM,method,side)
%
% Input:
%   x: signal to be forecasted
%   fs: sampling frequency
%   HOP: subsampling rate for forecasting
%   extK: size of the dataset for forecasting
%   extM: lengths of the segments used for forecasting
%   extSEC: extension length in seconds
%   method: forcasting method. To be chosen between
%       'SigExt': Least square estimation
%       'edmd': Empirical dynamical Mode decompostion
%       'gpr': Gaussian process regression
%       'symmetrization': Symmetric extension
%   side (optionnal): left blank for forward forecasting, set to 'backward' otherwise
% 
% Output:
%   xx: extended signal

L = round(fs*extSEC);

xext = forecasting(x,L,HOP,extK,extM,method);
xexti = forecasting(x,L,HOP,extK,extM,method,'backward');

xx = [xexti; x; xext];