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

p0 = 25 ;
phi00 = p0*fs/extM * ( tt + (0.3/(2*pi))*cos(2*pi*tt) ) ;
xx00 = cos(2*pi*phi00) ;

p1 = 60 ;
R = 1.4 + 0.2*cos(4*pi*tt) ;
phi01 = (p1*fs/extM) * tt + (250/2)*tt.^2 ;

nComp= 2;
switch nComp
    case 1
        xx01 = 0 ;
    case 2
        xx01 = R.*cos(2*pi*phi01) ;
end
xx0 = xx00 + xx01 ; % extended signal
x0 = xx0( (L+1) : (L+N) ) ; % restriction to the measurement interval

sigman = 1e-2 ;

%% Forecasting
nbXP = 1 ;

for ind = 1:nbXP
    noise = sigman*randn(N+2*L,1) ;
    x = x0.' + noise((L+1):(N+L)) ; % signal to be extended

    method.name = 'lse' ;
    tic;
    xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method).' ; % Forecasted signal via LSE
    LSEtime(ind) = toc;
    MeanLSE(ind,:) = xxLSE((N+L+1):end) - xx0((N+L+1):end) ;
    VarLSE(ind,:) = ( xxLSE((N+L+1):end) - xx0((N+L+1):end) ).^2 ;
    
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

BiasXP.LSE = mean(MeanLSE) ;
VarianceXP.LSE = mean(VarLSE) ;
CPUtimeXP.LSE = mean(LSEtime) ;

BiasXP.EDMD = mean(MeanEDMD) ;
VarianceXP.EDMD = mean(VarEDMD) ;
CPUtimeXP.EDMD = mean(EDMDtime) ;

BiasXP.GPR = mean(MeanGPR) ;
VarianceXP.GPR = mean(VarGPR) ;
CPUtimeXP.GPR = mean(GPRtime) ;


save('../../Results/PerfAHM','BiasXP','VarianceXP','CPUtimeXP');