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
        A = (Y*X') / (X*X') ; % least square estimation
        
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
        [Xi,mu,phix] = approxKoopman(X,Y,sigma2) ;
        
    case 'gpr'
        y = Y(end,:).' ;
        gprMdl = fitrgp(X.',y,'Basis','linear','FitMethod','exact','PredictMethod','exact');
        
end

% extension
L = round(fs*extSEC);
Z = zeros(extM,L) ;

switch method.name
    
    case 'edmd'
        tmp = phix.';
        for kk = 1:L
            tmp = mu.*tmp;
            Z(:,kk) =  tmp ; 
        end
        Z = real(Xi.' * Z);
        xext = Z(end,:)' ;
        
    case 'gpr'
        xext = zeros(L,1) ;
        z = x(end-extM+1: HOP: end).' ;
        for kk = 1:L
            zpred = predict(gprMdl,z);
            z = [z(2:end) zpred] ;
            xext(kk) = zpred ;
        end
        
    otherwise
        Z = zeros(extM,L) ;
        Z(:,1) = A*Y(:,end) ;
        for kk = 2:L
            Z(:,kk) = A*Z(:,kk-1) ; 
        end
        xext = Z(end,:)' ;
        
end

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    xext = flipud(xext);
end