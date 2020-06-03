clear all; 
%close all; clc;
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

tmp = cumsum(randn(N+2*L,1)) ;
tmp = smooth(tmp,0.1,'rloess') ;
tmp = tmp - mean(tmp) ;
tmp = tmp/max(abs(tmp));

p0 = 20 ;
freq00 = p0*fs/extM * (1 + 0.0*tmp) ;
% phi00 = p0*fs/extM * ( tt + (0.01/(2*pi))*cos(2*pi*tt) ) ;
phi00 = (1/fs) * cumsum(freq00).' ;
xx00 = cos(2*pi*phi00) ;


% phi01 = (p1*fs/extM) * tt + (20/2)*tt.^2 ;

nComp= 2;
switch nComp
    case 1
        xx01 = 0 ;
    case 2
        R = 1.4 + 0.2*cos(4*pi*tt) ;
        p1 = 35 ;
        tmp = cumsum(randn(N+2*L,1)) ;
        tmp = smooth(tmp,0.25,'rloess') ;
        tmp = tmp - mean(tmp) ;
        tmp = tmp/max(abs(tmp));
        freq01 = p1*fs/extM * (1 + 0.01*tmp) ;
        phi01 = (1/fs) * cumsum(freq01).' ;
        xx01 = R.*cos(2*pi*phi01) ;
end
xx0 = xx00 + xx01 ; % extended signal
x0 = xx0( (L+1) : (L+N) ) ; % restriction to the measurement interval

sigman = 1e-2 ;

%% Forecasting
nbXP = 1000 ;

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
MSE.LSE = BiasXP.LSE.^2 + VarianceXP.LSE ;
MSE.EDMD = BiasXP.EDMD.^2 + VarianceXP.EDMD ;
MSE.GPR = BiasXP.GPR.^2 + VarianceXP.GPR ;

fprintf('===================Forecasting Performance=========\n')
fprintf(' _________________________________________________\n')
fprintf('| Extension |  Computing  |   Forecasting Error   |\n')
fprintf('|  Method   | time (sec.) |    Mean   |  Variance |\n')
fprintf('|-----------|-------------|-----------|-----------|\n')
fprintf('|   LSE     |     %.2f    | %.3e | %.3e |\n', CPUtimeXP.LSE , median(MSE.LSE), var(MSE.LSE))
fprintf('|   EDMD    |     %.2f    | %.3e | %.3e |\n', CPUtimeXP.EDMD , median(MSE.EDMD), var(MSE.EDMD))
fprintf('|   GPR     |   %.2f    | %.3e | %.3e |\n', CPUtimeXP.GPR , median(MSE.GPR), var(MSE.GPR))
fprintf('|___________|_____________|___________|___________|\n\n')