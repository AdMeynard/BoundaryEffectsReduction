clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load ECG signal
addpath ../../Signals
data = load('PPG.mat');
x0 = data.PPG;

fs = data.fs ; % sampling frequency

N0 = length(x0) ;
N = 4e3; % subsignals length

% Forecasting parameters
methods = {'lse'} ;
HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

basicTF.hop = 10 ;
basicTF.win = 1501 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.01 ;

nXP = 12 ;
StdNoise = logspace(-3.5,-1,nXP) ;
for n = 1:nXP
    x = x0 + StdNoise(n)*randn(N0,1) ;
    TFR = {'sstSTFT','RS'} ;
    [forecastErr, CompTime, OTD] = SucForecast(x,fs,methods,HOP,N,extM,extK,extSEC,TFR,basicTF) ;
    
    NoiseForecastErr(n) = mean(forecastErr.LSE) ;
    NoiseTFR.SST(n) = mean(OTD.sst.LSE) ;
    NoiseTFR.STFT(n) = mean(OTD.stft.LSE) ;
    NoiseTFR.RS(n) = mean(OTD.rs.LSE) ;
end

SNR = 20*log10(std(x0)./StdNoise) ;

save('../../Results/resultSucForPPGnoise','SNR','NoiseForecastErr','NoiseTFR');