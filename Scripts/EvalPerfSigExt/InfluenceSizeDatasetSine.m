%% Evaluation, on a sine wave, of the influence of the size of the "training dataset" on SigExt's performance
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu


clear all; 
close all; clc;
addpath('../../Algorithm/');

%% Parameters

N = 10000 ; fs = N-1 ;
t = linspace(0,1,N);

% forecasting parameters
HOP = 1 ;
L = 100 ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length

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

%% Forecasting
sigman = 1e-2 ;

method.name = 'lse' ;
nbXP = 200 ;
nbXPP = 500 ;
KK = round(logspace(log10(2e2),log10(2e3),nbXP)) ;
Covh = zeros(extM,extM,L) ;
Ehw2 = zeros(L,1);
k = 1 ;
for k = 1:nbXP
    extK = KK(k) ;
    
    for nb2 = 1:nbXPP
        noise = sigman*randn(N+L,1) ;
        x = x0.' + noise(1:N) ; % signal to be extended
        
        % forecasting
        xext = forecasting(x,L,HOP,extK,extM,method).' ;
        BiasXP(nb2,:) = xext - xx0((N+1):end) ;
        MSEXP(nb2,:) = ( xext - xx0((N+1):end) ).^2 ;
    end
    BiasXPm(k,:) = mean(BiasXP) ;
    VarXPm(k,:) = mean(MSEXP) - BiasXPm(k,:).^2 ;
    
end

save('../../Results/PerfSizeDataset','BiasXPm','VarXPm','KK');