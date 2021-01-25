%% Evaluation of the Performance of SigExt on a AHM oscillatory signal
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; 
close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));

addpath('../../Algorithm/');

%% Parameters

N = 10000 ; fs = N-1 ;

% forecasting parameters
extSEC = 0.1 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;

t = linspace(0, (N-1)/fs, N) ;
tt = linspace(0, (N-1+L)/fs, N+L) ;

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

basicTF.fmin = 0 ;
basicTF.fmax = 0.05 ;
basicTF.df = (basicTF.fmax - basicTF.fmin)/1024 ;
basicTF.hop = 10 ; 
basicTF.win = 2*L+1 ;

tsstEXT = tt(1:basicTF.hop:end);

n0 = find(tsstEXT>=basicTF.win/2/fs,1,'first'); % avoid the left-side boundary
nf = find(tsstEXT<=max(t),1,'last');

for ind = 1:nbXP
    filename = strcat('../../Results/sigTBATSAHM/extTBATS_',num2str(ind),'.csv') ;
    xxTBATS = table2array( readtable(filename,'Range','B:B') ) ;
    
    x = xxTBATS(1:N) ;
    noise = sigman*randn(L,1) ;
    xx = [x; xx0((N+1):(N+L)).' + noise] ; % ground truth
    
    % TF reprensentations of the ground truth extension
    [~, SSTxx0, ~] = sqSTFT_C(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    [~, ~, conceFTxx0, ~] = ConceFT_sqSTFT_C(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    [STFTxx0, ~] = STFT_C(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    [~, RSxx0, ~, ~] = ConceFT_rsSTFT(xx, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    
    SSTxx0 = SSTxx0(:,n0:nf) ;
    conceFTxx0 = conceFTxx0(:,n0:nf) ;
    STFTxx0 = STFTxx0(:,n0:nf) ;
    RSxx0 = RSxx0(:,n0:nf) ;
    
    % TBATS

    
    [~, SSTxx, ~] = sqSTFT_C(xxTBATS, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    [~, ~, conceFTxx, ~] = ConceFT_sqSTFT_C(xxTBATS, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    [STFTxx, ~] = STFT_C(xxTBATS, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    [~, RSxx, ~, ~] = ConceFT_rsSTFT(xxTBATS, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    
    OTD.SST.TBATS(ind) = slicedOT(SSTxx(:,n0:nf), SSTxx0) ;
    OTD.conceFT.TBATS(ind) = slicedOT(conceFTxx(:,n0:nf), conceFTxx0) ;
    OTD.STFT.TBATS(ind) = slicedOT(STFTxx(:,n0:nf), STFTxx0) ;
    OTD.RS.TBATS(ind) = slicedOT(RSxx(:,n0:nf), RSxx0) ;
    
end

dataTBATS = table2array( readtable('../../Results/PerfAHM_TBATS.csv','Range','B:D','TreatAsEmpty','NA') ) ;
BiasXP.TBATS = dataTBATS(:,1) ;
VarianceXP.TBATS = dataTBATS(:,2) ;
CPUtimeXP.TBATS = dataTBATS(~isnan(dataTBATS(:,3)),3) ;


save('../../Results/PerfAHM_TBATS','BiasXP','VarianceXP','CPUtimeXP','OTD');
