%% Evaluation of the Performance of SigExt on a AHM oscillatory signal
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; 
close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));

addpath('../../Algorithm/');

%% Parameters

N = 10000 ; fs = N-1 ;
t = linspace(0,1,N);

% forecasting parameters
HOP = 1 ;
extSEC = 0.1 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = 750 ; % dimension of embedding / signals length
extK = round( 3.75*L );  % number of points to estimate A / size of datasets

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

sigman = 1e-2 ;

%% Forecasting + TFR
nbXP = 1000 ;

basicTF.fmin = 0 ;
basicTF.fmax = 0.05 ;
basicTF.df = (basicTF.fmax - basicTF.fmin)/1024 ;
basicTF.hop = 10 ; 
basicTF.win = 2*L+1 ;
LL = ceil(L/basicTF.hop); % extension in TF domain


t = linspace(0, (N-1)/fs, N) ;
tt = linspace(-L/fs, (N-1+L)/fs, N+2*L) ;
tsstEXT = tt(1:basicTF.hop:end);

n0 = find(tsstEXT>=min(t),1);
nf = find(tsstEXT<=max(t),1);

for ind = 1:nbXP
    noise = sigman*randn(N+2*L,1) ;
    xx = xx0.' + noise ; % ground truth
    
    x = xx((L+1):(N+L)) ; % signal to be extended
    
    [~, SSTxx0, ~] = sqSTFT_RT(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 0, 0) ;
    [~, ~, conceFTxx0, ~] = ConceFT_sqSTFT_RT(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1, 0, 0) ;
    [STFTxx0, ~] = STFT_RT(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10) ;
    [~, RSxx0, ~, ~] = ConceFT_rsSTFT_RT(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1) ;
    
    SSTxx0 = SSTxx0(:,n0:nf) ;
    conceFTxx0 = conceFTxx0(:,n0:nf) ;
    STFTxx0 = STFTxx0(:,n0:nf) ;
    RSxx0 = RSxx0(:,n0:nf) ;

    % SigExt
    method.name = 'SigExt' ;
    tic;
    xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; % Forecasted signal via SigExt
    LSEtime(ind) = toc;
    MeanLSE(ind,:) = xxLSE((N+L+1):end) - xx((N+L+1):end) ;
    VarLSE(ind,:) = ( xxLSE((N+L+1):end) - xx((N+L+1):end) ).^2 ;
    

    [~, SSTxx, ~] = sqSTFT_RT(xxLSE, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_RT(xxLSE, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_RT(xxLSE, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT_RT(xxLSE, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1) ;
    
    OTD.SST.LSE(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.LSE(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.LSE(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.LSE(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
    % Symmetrization
    method.name = 'symmetrization' ;
    tic;
    xxSYM = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ;
    SYMtime(ind) = toc;
    MeanSYM(ind,:) = xxSYM((N+L+1):end) - xx((N+L+1):end) ;
    VarSYM(ind,:) = ( xxSYM((N+L+1):end) - xx((N+L+1):end) ).^2 ;
    
    [~, SSTxx, ~] = sqSTFT_RT(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_RT(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_RT(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT_RT(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1) ;
    
    OTD.SST.SYM(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.SYM(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.SYM(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.SYM(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
    % EDMD
    method.name = 'edmd' ;
    method.param = 100 ;
    tic; 
    xxEDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; 
    EDMDtime(ind) = toc;
    MeanEDMD(ind,:) = xxEDMD((N+L+1):end) - xx((N+L+1):end) ;
    VarEDMD(ind,:) = ( xxEDMD((N+L+1):end) - xx((N+L+1):end) ).^2 ;
    
    [~, SSTxx, ~] = sqSTFT_RT(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_RT(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_RT(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT_RT(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1) ;
    
    OTD.SST.EDMD(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.EDMD(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.EDMD(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.EDMD(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
    % GPR
    method.name = 'gpr' ;
    tic; 
    xxGPR = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; 
    GPRtime(ind) = toc;
    MeanGPR(ind,:) = xxGPR((N+L+1):end) - xx((N+L+1):end) ;
    VarGPR(ind,:) = ( xxGPR((N+L+1):end) - xx((N+L+1):end) ).^2 ;
    
    [~, SSTxx, ~] = sqSTFT_RT(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_RT(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_RT(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT_RT(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1) ;
    
    OTD.SST.GPR(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.GPR(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.GPR(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.GPR(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
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

save('../../Results/PerfAHM','BiasXP','VarianceXP','CPUtimeXP','OTD');
