function xext = forecasting(x,fs,HOP,extN,extP,extSEC,varargin)

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    x = flipud(x);
end

X = [] ; Y = [] ;
for kk = 1: extN
    X(:,kk) = x(end-extN-extP+kk: HOP: end-extN+kk-1) ;
    Y(:,kk) = x(end-extN-extP+kk+1: HOP: end-extN+kk) ;
end


% A is extP * extP
A = Y*X'*inv(X*X') ; % least square estimation

% extension
Z = [] ;
Z(:,1) = Y(:,end) ;

for kk = 1:fs*extSEC
    Z(:,kk+1) = A*Z(:,kk) ; 
end
xext = Z(end, 2:end)' ;

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    xext = flipud(xext);
end