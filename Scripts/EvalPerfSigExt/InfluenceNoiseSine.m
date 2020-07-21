%% Evaluation of Gaussian white noise on SigExt's performance on a sine wave
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; 
close all; clc;
addpath('../../Algorithm/');

%% Parameters

% forecasting parameters
HOP = 1 ;
L = 350 ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 3*extM );  % number of points to estimate A / size of datasets

N = extK + extM + 1 ; 
fs = N-1 ;
t = linspace(0,1,N);
tt = linspace(0, 1+L/fs, N+L) ;

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
x0 = xx0(1:N) ; % restriction to the measurement interval

%% Ideal Forecasting vector

X0 = [] ; Y0 = [] ;
for kk = 1: extK
    X0(:,kk) = x0(end-extK-extM+kk: HOP: end-extK+kk-1) ;
    Y0(:,kk) = x0(end-extK-extM+kk+1: HOP: end-extK+kk) ;
end
zK = x0(end-extM+1 : HOP: end).' ;

%% Forecasting
method.name = 'lseV' ;
nbXP = 3 ; % number of noise levels
nbXPP = 3000 ; % number of XP / noise level
Sigma = logspace(-3,-1,nbXP) ;

hopL = 10 ; % hop for forecasting
lenForcast = length(1:hopL:L) ;
Covh = zeros(extM,extM, lenForcast) ;
Ehw2 = zeros(lenForcast,1);
k = 1 ;
for k = 1:nbXP
    sigman = Sigma(k) ;
    
    % Theoretical foreacstion Matrix
    S0 = (1/extK)*(X0*X0') + sigman^2*eye(extM) ;
    r = zeros(1,extM) ; r(2) = 1 ;
    S1 = (1/extK)*(Y0*X0') + sigman^2*toeplitz(zeros(extM,1),r) ;
    A0 = S1 / S0 ; % ideal least square estimation
    j = 1;
    for l = 1:hopL:L
        A0l{j} = A0^l ;
        alpha0l(j,:) = A0l{j}(end,:) ;
        j = j + 1 ;
    end
    a2 = sum(alpha0l.^2, 2) ;
    
    for nb2 = 1:nbXPP
        noise = sigman*randn(N+L,1) ;
        wK = noise((N-extM+1):N)/sigman ;
        x = x0.' + noise(1:N) ; % signal to be extended
        
        % Evaluate h
        X = [] ; Y = [] ;
        for kk = 1: extK
            X(:,kk) = x(end-extK-extM+kk: HOP: end-extK+kk-1) ;
            Y(:,kk) = x(end-extK-extM+kk+1: HOP: end-extK+kk) ;
        end
        A = (Y*X') / (X*X') ; % experimental least square estimation
        j = 1 ;
        for l = 1:hopL:L
            tmp = A^l - A0l{j} ;
            h = tmp(end,:) ;
            Covh(:,:,j) = Covh(:,:,j) + h' * h ;
            Ehw2(j) = Ehw2(j) + (h*wK)^2 ; % !!! TESTS
            j  = j + 1 ;
        end
        
        % forecasting
        xext = forecasting(x,L,HOP,extK,extM,method).' ; %SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ;
        MeanXP(nb2,:) = xext - xx0((N+1):end) ;
        VarXP(nb2,:) = ( xext - xx0((N+1):end) ).^2 ;
    end
    MeanXPm(k,:) = mean(MeanXP) ;
    VarXPm(k,:) = mean(VarXP) ;
    
    Gammal = extK*Covh/nbXPP ;
    Ehw2 = Ehw2/nbXPP ;
    j = 1 ;
    for l=1:hopL:L
        VarTH(k,j) = (1/extK) * zK'*Gammal(:,:,j)*zK + sigman^2*a2(j) +  sigman^2*Ehw2(j); %(sigman^2/extK) * trace(Gammal(:,:,l)) ;
        j  = j +  1;
    end
end

save('../../Results/PerfNoise','MeanXPm','VarXPm','VarTH','extM','Sigma');