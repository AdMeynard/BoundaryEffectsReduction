clear all; close all; clc;
addpath(genpath('../../SST'));
addpath('../../algorithm/');

%% load ECG signal
addpath ../../Signals
data = load('ECG.mat');
xtot = data.ECG;

fs = 200 ; % sampling frequency

N = 15e3; % subsignals length

% Forecasting parameters
HOP = 1 ;
extSEC = 6 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

basicTF.meth = 'sstSTFT' ;
basicTF.hop = 10 ;
basicTF.win = 2301 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.01 ;
[forecastErr1, CompTime1, OTD1] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,basicTF) ;
basicTF.meth = 'RS' ;
[forecastErr2, CompTime2, OTD2] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,basicTF) ;


save('../../Results/resultSucForECGcompMeth','forecastErr1','CompTime1','OTD1','forecastErr2','CompTime2','OTD2');