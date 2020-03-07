clear all; 
%close all; clc;
addpath('../../algorithm/');

%% Parameters

N = 10000 ; fs = N-1 ;
t = linspace(0,1,N);

f = 100;
a = 0;%100;
x0 = cos(2*pi*(0.5*a*t.^2+f*t));

% forecasting parameters
HOP = 1 ;
extSEC = 0.05 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets
if extK + extM >length(x0) - 10 
    extK = extK/2 ; extM = extM/2 ;
end

tt = linspace(-L/fs, 1+L/fs, N+2*L) ;
xx0L = cos(2*pi*(0.5*a*tt.^2+f*tt));

%% Forecasting
method.name = 'lseV' ;
nbXP = 100 ;
nbXPP = 500;
Sigma = linspace(5e-3,1e-1,nbXP) ;
MSE10m = zeros(nbXP,1) ;
MSE100m = zeros(nbXP,1) ;
MSE200m = zeros(nbXP,1) ;
MSE500m = zeros(nbXP,1) ;
k = 1 ;
for sigman = Sigma
    for nb2 = 1:nbXPP
        noise = sigman*randn(N+2*L,1) ;
        x = x0.' + noise((L+1):(N+L)) ;
        xxTRUE = xx0L.' + noise ;
        
        xx = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ;
        MSE10(nb2) = (xx(N+L+10) - xxTRUE(N+L+10))^2 ;
        MSE100(nb2) = (xx(N+L+100) - xxTRUE(N+L+100))^2 ;
        MSE200(nb2) = (xx(N+L+200) - xxTRUE(N+L+200))^2 ;
        MSE500(nb2) = (xx(N+L+500) - xxTRUE(N+L+500))^2 ;
    end
    MSE10m(k) = mean(MSE10) ;
    MSE100m(k) = mean(MSE100) ;
    MSE200m(k) = mean(MSE200) ;
    MSE500m(k) = mean(MSE500) ;
    k = k+1 ;
end

save('../../Results/MSEnoise','MSE10m','MSE100m','MSE200m','MSE500m');