%% Evaluate BoundEffRed on a EEG signal
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load ECG signal
addpath ../../Signals
data = load('EEG.mat');
xtot = data.EEG;

fs = 200 ; % sampling frequency

N = 3e3; % subsignals length

% Forecasting parameters
methods = {'SigExt','symmetrization','edmd','gpr'} ;
HOP = 1 ;
extSEC = 1.5 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

TFR = {'conceFT','RS'} ;
basicTF.hop = 10 ;
basicTF.win = 2*L+1 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.15 ;
basicTF.df = 2e-5 ;
basicTF.MT = 30 ;
[forecastErr, CompTime, OTD] = SucForecast(xtot,fs,methods,HOP,N,extM,extK,extSEC,TFR,basicTF) ;

save('../../Results/resultSucForEEGcompMeth','forecastErr','CompTime','OTD');
