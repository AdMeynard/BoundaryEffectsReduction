clear all; close all; clc;
addpath(genpath('../../SST'));
addpath('../../algorithm/');

%% load ECG signal
addpath ../../Signals
data = load('PPGsig.mat');
xtot = data.PPG;

fs = data.fs ; % sampling frequency

N = 4e3; % subsignals length

% Forecasting parameters
HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

basicTF.meth = 'sstSTFT' ;
basicTF.hop = 10 ;
basicTF.win = 1501 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.01 ;
[forecastErr1, sstS, sstLSE, stftS, stftLSE] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,basicTF) ;
basicTF.meth = 'RS' ;
[forecastErr2, rsS, rsLSE] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,basicTF) ;


save('../../Results/resultSucForPPGcompMeth','forecastErr1','sstS','sstLSE','stftS', 'stftLSE', 'forecastErr2', 'rsS', 'rsLSE');