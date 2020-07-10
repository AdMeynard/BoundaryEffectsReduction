clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load ECG signal
addpath ../../Signals
data = load('PPG.mat');
xtot = data.PPG;

fs = data.fs ; % sampling frequency

N = 4e3; % subsignals length

% Forecasting parameters
methods = {'lse','symmetrization','edmd','gpr'} ;
HOP = 1 ;
extSEC = 5.2 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

TFR = {'sstSTFT','conceFT'} ;
basicTF.hop = 10 ;
basicTF.win = 1501 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.01 ;
[forecastErr, CompTime, OTD] = SucForecast(xtot,fs,methods,HOP,N,extM,extK,extSEC,TFR,basicTF) ;

save('../../Results/resultSucForPPGcompMeth','forecastErr','CompTime','OTD');