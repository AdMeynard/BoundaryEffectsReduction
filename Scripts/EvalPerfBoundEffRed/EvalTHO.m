clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load ECG signal
addpath ../../Signals
data = load('THO.mat');
xtot = data.THO;

fs = 100 ; % sampling frequency

N = 6e3; % subsignals length

% Forecasting parameters
HOP = 1 ;
extSEC = 7 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

basicTF.meth = 'sstSTFT' ;
basicTF.hop = 10 ;
basicTF.win = 1501 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.01 ;
[forecastErr1, CompTime1, OTD1] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,basicTF) ;
basicTF.meth = 'RS' ;
[forecastErr2, CompTime2, OTD2] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,basicTF) ;


save('../../Results/resultSucForTHOcompMeth','forecastErr1','CompTime1','OTD1','forecastErr2','CompTime2','OTD2');