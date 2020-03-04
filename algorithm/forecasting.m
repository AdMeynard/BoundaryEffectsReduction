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
switch method.name
    
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
        
    case 'edmd'
        sigma2 = method.param ;
        [Xi,mu,phix] = approxKoopmanKB(X,Y,sigma2) ;
        
end

% extension
K = round(fs*extSEC);
Z = zeros(extM,K) ;

switch method.name
    
    case 'edmd'
        tmp = phix.';
        for kk = 1:K
            tmp = mu.*tmp;
            Z(:,kk) =  tmp ; 
        end
        Z = real(Xi.' * Z);
        
    otherwise
        Z = zeros(extM,K) ;
        Z(:,1) = A*Y(:,end) ;
        for kk = 2:K
            Z(:,kk) = A*Z(:,kk-1) ; 
        end
        
end

xext = Z(end,:)' ;

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    xext = flipud(xext);
end