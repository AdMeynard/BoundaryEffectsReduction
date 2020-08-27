%% Evaluation of Gaussian white noise on SigExt's performance on a sine wave
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; 
close all; clc; warning off;
addpath('../../Algorithm/');

%% Parameters

% forecasting parameters
HOP = 1 ;
L = 100 ;
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
method.name = 'SigExt' ;
nbXP = 200 ; % number of noise levels
nbXPP = 500 ; % number of XP / noise level
Sigma = logspace(-7,-1,nbXP) ;

hopL = 10 ; % hop for forecasting
lenForcast = length(1:hopL:L) ;
Covh = zeros(extM,extM, lenForcast) ;
Ehw2 = zeros(lenForcast,1);
k = 1 ;
for k = 1:nbXP
    sigman = Sigma(k) ;
    
    for nb2 = 1:nbXPP
        noise = sigman*randn(N+L,1) ;
        x = x0.' + noise(1:N) ; % signal to be extended
        
        % forecasting
        xext = forecasting(x,L,HOP,extK,extM,method).' ;
        BiasXP(nb2,:) = xext - xx0((N+1):end) ;
        MSE_XP(nb2,:) = ( xext - xx0((N+1):end) ).^2 ;
    end
    
    BiasXPm(k,:) = mean(BiasXP) ;
    VarXPm(k,:) = mean(MSE_XP) - BiasXPm(k,:).^2 ;
    
end

save('../../Results/PerfNoise','BiasXPm','VarXPm','Sigma');