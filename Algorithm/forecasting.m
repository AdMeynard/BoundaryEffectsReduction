function xext = forecasting(x,L,HOP,extK,extM,method,varargin)
% FORECASTING Predict future values of a signal
% Usage:	xext = forecasting(x,L,HOP,extK,extM,method,side)
%
% Input:
%   x: signal to be forecasted
%   L: number of samples to be forecasted
%   HOP: subsampling rate for forecasting
%   extK: size of the dataset for forecasting
%   extM: lengths of the segments used for forecasting
%   method: forcasting method. To be chosen between
%       'SigExt': Least square estimation
%       'edmd': Empirical dynamical Mode decompostion
%       'gpr': Gaussian process regression
%       'symmetrization': Symmetric extension
%   side (optionnal): left blank for forward forecasting, set to 'backward' otherwise
% 
% Output:
%   xext: forecasted extension

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    x = flipud(x);
end

if ~strcmpi(method.name,'symmetrization')
    X = zeros(extM,extK) ;
    for kk = 1: extK
        X(:,kk) = x(end-extK-extM+kk: HOP: end-extK+kk-1) ;
    end
    Y = [X(:,2:end) x(end-extM+1: HOP: end)] ;
end


%% Estimate the parameters of the forecasting model
switch method.name
    
    case {'SigExt','lse'}
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
        
    otherwise
        % nothing
end

%% Extension
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
        
    case {'SigExt','lse','lseV'}
        Z = zeros(extM,L) ;
        Z(:,1) = A*Y(:,end) ;
        for kk = 2:L
            Z(:,kk) = A*Z(:,kk-1) ; 
        end
        xext = Z(end,:)' ;
        
    otherwise
        xext = flipud( x((end-L+1):end) ) ; % symmetrization
        
end

if (~isempty(varargin))&&(isequal(varargin{1},'backward'))
    xext = flipud(xext);
end