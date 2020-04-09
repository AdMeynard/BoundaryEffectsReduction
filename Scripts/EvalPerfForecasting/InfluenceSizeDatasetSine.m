clear all; 
%close all; clc;
addpath('../../Algorithm/');

%% Parameters

N = 10000 ; fs = N-1 ;
t = linspace(0,1,N);

% forecasting parameters
HOP = 1 ;
extSEC = 0.05 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets

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

sigman = 1e-2 ;

%% Forecasting
method.name = 'lseV' ;
nbXP = 50 ;
nbXPP = 3000 ;
KK = round(linspace(1e3,6e3,nbXP)) ;
k = 1 ;
for k = 1:nbXP
    extK = KK(k) ;
    for nb2 = 1:nbXPP
        noise = sigman*randn(N+2*L,1) ;
        x = x0.' + noise((L+1):(N+L)) ; % signal to be extended
        
        xx = SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ;
        MeanXP(nb2,:) = xx((N+L+1):end) - xx0((N+L+1):end) ;
        VarXP(nb2,:) = ( xx((N+L+1):end) - xx0((N+L+1):end) ).^2 ;
    end
    MeanXPm(k,:) = mean(MeanXP) ;
    VarXPm(k,:) = mean(VarXP) ;
end

save('../../Results/PerfSizeDataset','MeanXPm','VarXPm','extM','Sigma');