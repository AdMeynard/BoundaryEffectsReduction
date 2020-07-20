%% Evaluation the "training dataset" on SigExt's performance on a sine wave
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu


clear all; 
%close all; clc;
addpath('../../Algorithm/');

%% Parameters

N = 10000 ; fs = N-1 ;
t = linspace(0,1,N);

% forecasting parameters
HOP = 1 ;
extSEC = 0.01 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length

tt = linspace(-L/fs, 1+L/fs, N+2*L) ;

%% Synthesize signal

p0 = 10 ;
f0 = p0*fs/extM ;
xx00 = cos(2*pi*f0*tt) ;

p1 = 33 ;
R = 1.4 ;
f1 = p1*fs/extM ;

nComp= 2;
switch nComp
    case 1
        xx01 = 0 ;
    case 2
        xx01 = R*cos(2*pi*f1*tt) ;
end
xx0 = xx00 + xx01 ; % extended signal
x0 = xx0( (L+1) : (L+N) ) ; % restriction to the measurement interval

zK = x0(end-extM+1 : HOP: end).' ;

%% Forecasting
sigman = 1e-2 ;

method.name = 'lse' ;
nbXP = 3 ;
nbXPP = 3000 ;
KK = round(logspace(log10(4.5e2),log10(2e3),nbXP)) ;
Covh = zeros(extM,extM,L) ;
Ehw2 = zeros(L,1);
k = 1 ;
for k = 1:nbXP
    extK = KK(k) ;
    
    X0 = [] ; Y0 = [] ;
    for kk = 1: extK
        X0(:,kk) = x0(end-extK-extM+kk: HOP: end-extK+kk-1) ;
        Y0(:,kk) = x0(end-extK-extM+kk+1: HOP: end-extK+kk) ;
    end
    S0 = (1/extK)*(X0*X0') + sigman^2*eye(extM) ;
    r = zeros(1,extM) ; r(2) = 1 ;
    S1 = (1/extK)*(Y0*X0') + sigman^2*toeplitz(zeros(extM,1),r) ;
    A0 = S1 / S0 ; % ideal least square estimation
    for l = 1:L
        tmp0 = A0^l ;
        alpha0l(l,:) = tmp0(end,:) ;
    end
    a2 = sum(alpha0l.^2, 2) ;
    
    for nb2 = 1:nbXPP
        noise = sigman*randn(N+2*L,1) ;
        wK = noise((N+L-extM+1):(N+L))/sigman ;
        x = x0.' + noise((L+1):(N+L)) ; % signal to be extended
        
                % Evaluate h
        X = [] ; Y = [] ;
        for kk = 1: extK
            X(:,kk) = x(end-extK-extM+kk: HOP: end-extK+kk-1) ;
            Y(:,kk) = x(end-extK-extM+kk+1: HOP: end-extK+kk) ;
        end
        A = (Y*X') / (X*X') ; % experimental least square estimation
        for l = 1:L
            tmp = A^l - A0^l ;
            h = tmp(end,:) ;
            Covh(:,:,l) = Covh(:,:,l) + h' * h ;
            Ehw2(l) = Ehw2(l) + (h*wK)^2 ; % !!! TESTS
        end
        
        % forecasting
        xx = SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ;
        MeanXP(nb2,:) = xx((N+L+1):end) - xx0((N+L+1):end) ;
        VarXP(nb2,:) = ( xx((N+L+1):end) - xx0((N+L+1):end) ).^2 ;
    end
    MeanXPm(k,:) = mean(MeanXP) ;
    VarXPm(k,:) = mean(VarXP) ;
    
    Gammal = extK*Covh/nbXPP ;
    Ehw2 = Ehw2/nbXPP ;
    for l=1:L
        VarTH(k,l) = (1/extK) * zK'*Gammal(:,:,l)*zK + sigman^2*a2(l) +  sigman^2*Ehw2(l); %(sigman^2/extK) * trace(Gammal(:,:,l)) ;
    end
end

save('../../Results/PerfSizeDataset','MeanXPm','VarXPm','VarTH','extM','sigman','KK');