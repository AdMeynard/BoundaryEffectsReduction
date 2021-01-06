%% Evaluation on a AHM signal of the sunsignals length influence on SigExt's performance
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
extM = 750 ;
extK = round( 3.75*L );  % number of points to estimate A / size of datasets

tt = linspace(0, 1+L/fs, N+L) ;

%% Synthesize signal

M0 = 750 ;
p0 = 10 ;
phi00 = p0/M0 *  ( tt*fs + (0.01/(2*pi))*cos(2*pi*tt*fs/N) ) ;
xx00 = cos(2*pi*phi00) ;

nComp= 2;
switch nComp
    case 1
        xx01 = 0 ;
    case 2
        R = 1.4 + 0.2*cos(4*pi*tt) ;
        p1 = 23 ;
        phi01 = (p1*fs/M0) * tt + (20/2)*tt.^2 ;
        xx01 = R.*cos(2*pi*phi01) ;
end
xx0 = xx00 + xx01 ; % extended signal

sigman = 1e-2 ;

%% Forecasting
nbXP = 1000 ;


for indXP = 1:nbXP
    noise = sigman*randn(N+L,1) ;
    xx = xx0.' + noise ;
    
    x = xx(1:N) ; % signal to be extended

    method.name = 'symmetrization' ;
    tic;
    xxSYM = forecasting(x,L,HOP,extK,extM,method) ; % Forecasted signal via SigExt  
    CPUtimeXP.SYM(indXP) = toc;
    Bias.SYM(indXP) = mean( xxSYM - xx((N+1):end) ) ;
    MSE.SYM(indXP) = mean( ( xxSYM - xx((N+1):end) ).^2 );
    
    method.name = 'edmd' ;
    method.param = 200 ;
    tic;
    xxEDMD = forecasting(x,L,HOP,extK,extM,method) ; % Forecasted signal via SigExt  
    CPUtimeXP.EDMD(indXP) = toc;
    Bias.EDMD(indXP) = mean( xxEDMD - xx((N+1):end) ) ;
    MSE.EDMD(indXP) = mean( ( xxEDMD - xx((N+1):end) ).^2 ) ;
    
    method.name = 'gpr' ;
    tic;
    xxGPR = forecasting(x,L,HOP,extK,extM,method) ; % Forecasted signal via SigExt  
    CPUtimeXP.GPR(indXP) = toc;
    Bias.GPR(indXP) = mean( xxGPR - xx((N+1):end) ) ;
    MSE.GPR(indXP) = mean( ( xxGPR - xx((N+1):end) ).^2 ) ;
    
end

save('../../Results/PerfOtherMethodsAHM','Bias','MSE','CPUtimeXP');

