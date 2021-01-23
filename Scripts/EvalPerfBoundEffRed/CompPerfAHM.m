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
extK = round( 3.75*L ) ;  % number of points to estimate A / size of datasets

tt = linspace(0, 1+L/fs, N+L) ;

%% Synthesize signal

p0 = 10 ;
M0 = 750 ;
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

%% Forecasting + TFR
nbXP = 100 ;

extMval = [100 750 1500 3000]; % dimension of embedding / signals length

basicTF.fmin = 0 ;
basicTF.fmax = 0.05 ;
basicTF.df = (basicTF.fmax - basicTF.fmin)/1024 ;
basicTF.hop = 10 ; 
basicTF.win = 2*L+1 ;


t = linspace(0, (N-1)/fs, N) ;
tt = linspace(0, (N-1+L)/fs, N+L) ;
tsstEXT = tt(1:basicTF.hop:end);

n0 = find(tsstEXT>=basicTF.win/2/fs,1,'first'); % avoid the left-side boundary
nf = find(tsstEXT<=max(t),1,'last');

for ind = 1:nbXP
    noise = sigman*randn(N+L,1) ;
    xx = xx0.' + noise ; % ground truth
    
    % TF reprensentations of the ground truth extension
    [~, SSTxx0, ~] = sqSTFT_C(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    [~, ~, conceFTxx0, ~] = ConceFT_sqSTFT_C(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    [STFTxx0, ~] = STFT_C(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    [~, RSxx0, ~, ~] = ConceFT_rsSTFT(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    
    SSTxx0 = SSTxx0(:,n0:nf) ;
    conceFTxx0 = conceFTxx0(:,n0:nf) ;
    STFTxx0 = STFTxx0(:,n0:nf) ;
    RSxx0 = RSxx0(:,n0:nf) ;
    
    % Signal to be extended
    x = xx(1:N) ; 
    [~, SSTx, ~] = sqSTFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    [~, ~, conceFTx, ~] = ConceFT_sqSTFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    [STFTx, ~] = STFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    [~, RSx, ~, ~] = ConceFT_rsSTFT(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    
    OTD.SST.S(ind) = slicedOT(SSTx(:,n0:end), SSTxx0) ;
    OTD.conceFT.S(ind) = slicedOT(conceFTx(:,n0:end), conceFTxx0) ;
    OTD.STFT.S(ind) = slicedOT(STFTx(:,n0:end), STFTxx0) ;
    OTD.RS.S(ind) = slicedOT(RSx(:,n0:end), RSxx0) ;

    % SigExt
    method.name = 'SigExt' ;
    indM = 1 ;
    for extM = extMval
        tic;
        xxSigExt = forecasting(x,L,HOP,extK,extM,method) ; % Forecasted signal via SigExt  
        LSEtime(ind,indM) = toc;
        MeanLSE(ind,indM) = mean( xxSigExt - xx((N+1):end) ) ;
        VarLSE(ind,indM) = mean( ( xxSigExt - xx((N+1):end) ).^2 );
        
        xxSigExt = [x; xxSigExt] ;
        [~, SSTxx, ~] = sqSTFT_C(xxSigExt, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
        [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_C(xxSigExt, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        [STFTxx, ~] = STFT_C(xxSigExt, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
        [~, RSxx, ~, ~] = ConceFT_rsSTFT(xxSigExt, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;

        OTD.SST.LSE(ind,indM) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
        OTD.conceFT.LSE(ind,indM) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
        OTD.STFT.LSE(ind,indM) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
        OTD.RS.LSE(ind,indM) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
        
        indM = indM+1 ;
    end
    
    % Symmetrization
    method.name = 'symmetrization' ;
    tic;
    xxSYM = forecasting(x,L,HOP,extK,extM,method) ;
    SYMtime(ind) = toc;
    MeanSYM(ind,:) = xxSYM - xx((N+1):end) ;
    VarSYM(ind,:) = ( xxSYM - xx((N+1):end) ).^2 ;
    
    xxSYM = [x; xxSYM] ;
    [~, SSTxx, ~] = sqSTFT_C(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_C(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_C(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT(xxSYM, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    
    OTD.SST.SYM(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.SYM(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.SYM(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.SYM(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
    % EDMD
    method.name = 'edmd' ;
    method.param = 100 ;
    tic; 
    xxEDMD = forecasting(x,L,HOP,extK,extM,method) ;
    EDMDtime(ind) = toc;
    MeanEDMD(ind,:) = xxEDMD - xx((N+1):end) ;
    VarEDMD(ind,:) = ( xxEDMD - xx((N+1):end) ).^2 ;
    
    xxEDMD = [x; xxEDMD] ;
    [~, SSTxx, ~] = sqSTFT_C(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_C(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_C(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT(xxEDMD, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    
    OTD.SST.EDMD(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.EDMD(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.EDMD(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.EDMD(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
    % GPR
    method.name = 'gpr' ;
    tic; 
    xxGPR = forecasting(x,L,HOP,extK,extM,method) ;
    GPRtime(ind) = toc;
    MeanGPR(ind,:) = xxGPR - xx((N+1):end) ;
    VarGPR(ind,:) = ( xxGPR - xx((N+1):end) ).^2 ;
    
    xxGPR = [x; xxGPR] ;
    [~, SSTxx, ~] = sqSTFT_C(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_C(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_C(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT(xxGPR, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    
    OTD.SST.GPR(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.GPR(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.GPR(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.GPR(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
end

BiasXP.LSE = MeanLSE ;
VarianceXP.LSE = VarLSE ;
CPUtimeXP.LSE = mean(LSEtime,1) ;

BiasXP.SYM = mean(MeanSYM,2) ;
VarianceXP.SYM = mean(VarSYM,2) ;
CPUtimeXP.SYM = mean(SYMtime) ;

BiasXP.EDMD = mean(MeanEDMD,2) ;
VarianceXP.EDMD = mean(VarEDMD,2) ;
CPUtimeXP.EDMD = mean(EDMDtime) ;

BiasXP.GPR = mean(MeanGPR,2) ;
VarianceXP.GPR = mean(VarGPR,2) ;
CPUtimeXP.GPR = mean(GPRtime) ;

save('../../Results/PerfAHM','BiasXP','VarianceXP','CPUtimeXP','OTD','extMval');
