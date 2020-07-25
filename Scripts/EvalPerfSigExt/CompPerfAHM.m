%% Evaluation of the Performance of SigExt on a AHM oscillatory signal
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
extSEC = 0.1 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets

tt = linspace(-L/fs, 1+L/fs, N+2*L) ;

%% Synthesize signal

p0 = 10 ;
phi00 = p0/extM *  ( tt*fs + (0.01/(2*pi))*cos(2*pi*tt*fs/N) ) ;
xx00 = cos(2*pi*phi00) ;

nComp= 2;
switch nComp
    case 1
        xx01 = 0 ;
    case 2
        R = 1.4 + 0.2*cos(4*pi*tt) ;
        p1 = 23 ;
        phi01 = (p1*fs/extM) * tt + (20/2)*tt.^2 ;
        xx01 = R.*cos(2*pi*phi01) ;
end
xx0 = xx00 + xx01 ; % extended signal
x0 = xx0( (L+1) : (L+N) ) ; % restriction to the measurement interval

sigman = 1e-2 ;

%% Forecasting
nbXP = 50 ;

for ind = 1:nbXP
    noise = sigman*randn(N+2*L,1) ;
    x = x0.' + noise((L+1):(N+L)) ; % signal to be extended

    method.name = 'SigExt' ;
    tic;
    xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ; % Forecasted signal via SigExt
    LSEtime(ind) = toc;
    MeanLSE(ind,:) = xxLSE((N+L+1):end) - xx0((N+L+1):end) ;
    VarLSE(ind,:) = ( xxLSE((N+L+1):end) - xx0((N+L+1):end) ).^2 ;
    
    method.name = 'symmetrization' ;
    tic;
    xxSYM = SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ;
    SYMtime(ind) = toc;
    MeanSYM(ind,:) = xxSYM((N+L+1):end) - xx0((N+L+1):end) ;
    VarSYM(ind,:) = ( xxSYM((N+L+1):end) - xx0((N+L+1):end) ).^2 ;
    
    method.name = 'edmd' ;
    method.param = 100 ;
    tic; 
    xxEDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ; 
    EDMDtime(ind) = toc;
    MeanEDMD(ind,:) = xxEDMD((N+L+1):end) - xx0((N+L+1):end) ;
    VarEDMD(ind,:) = ( xxEDMD((N+L+1):end) - xx0((N+L+1):end) ).^2 ;
    
    method.name = 'gpr' ;
    tic; 
    xxGPR = SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ; 
    GPRtime(ind) = toc;
    MeanGPR(ind,:) = xxGPR((N+L+1):end) - xx0((N+L+1):end) ;
    VarGPR(ind,:) = ( xxGPR((N+L+1):end) - xx0((N+L+1):end) ).^2 ;
    
end

BiasXP.LSE = mean(MeanLSE,2) ;
VarianceXP.LSE = mean(VarLSE,2) ;
CPUtimeXP.LSE = mean(LSEtime) ;

BiasXP.SYM = mean(MeanSYM,2) ;
VarianceXP.SYM = mean(VarSYM,2) ;
CPUtimeXP.SYM = mean(SYMtime) ;

BiasXP.EDMD = mean(MeanEDMD,2) ;
VarianceXP.EDMD = mean(VarEDMD,2) ;
CPUtimeXP.EDMD = mean(EDMDtime) ;

BiasXP.GPR = mean(MeanGPR,2) ;
VarianceXP.GPR = mean(VarGPR,2) ;
CPUtimeXP.GPR = mean(GPRtime) ;

save('../../Results/PerfAHM','BiasXP','VarianceXP','CPUtimeXP');
