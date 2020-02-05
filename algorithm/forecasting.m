function xext = forecasting(x,fs,HOP,extK,extM,extSEC,method,varargin)

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    x = flipud(x);
end

X = [] ; Y = [] ;
for kk = 1: extK
    X(:,kk) = x(end-extK-extM+kk: HOP: end-extK+kk-1) ;
    Y(:,kk) = x(end-extK-extM+kk+1: HOP: end-extK+kk) ;
end


% A is extP * extP
switch method
    
    case 'lse'
        A = Y*X'*inv(X*X') ; % least square estimation
        
    case 'lseV'
        r = zeros(1,extM); r(2)=1;
        B = toeplitz(zeros(extM,1),r);
        eMeMt = zeros(extM,extM); eMeMt(extM,extM)=1;
        A = B + (eMeMt*Y*X')/(X*X') ; % vector least square estimation

    case 'dmd'
        [U,S,V] = svd(X,'econ');
        iS = inv(S);
        A = Y*V*iS*U';
        
end

% extension
Z = [] ;
Z(:,1) = Y(:,end) ;

K = round(fs*extSEC);
for kk = 1:K
    Z(:,kk+1) = A*Z(:,kk) ; 
end
xext = Z(end, 2:end)' ;

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    xext = flipud(xext);
end