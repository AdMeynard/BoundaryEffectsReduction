%% Real-time SST using BoundEffRed
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');
addpath('../../Signals/');

%% load THO signal
data = load('PPG.mat');

subT = 2 ; % Sub sampling (pre-processing)
xtot = data.PPG(1:subT:end);
fs = data.fs/subT ; % sampling frequency

minutes = 3 ; 
Nmax = minutes*60*fs ;
xtot = xtot(1:Nmax) ; % shorten signal

%% Parameters
forecastMethod.name = 'SigExt' ;

basicTF.representation = 'SST' ;

basicTF.hop = round(10/subT) ;
basicTF.win = 2*round(1000/2/subT)+1 ; % window length (in samples) /!\ must be odd
basicTF.fmin = 0.5/fs ;
basicTF.fmax = 5.5/fs ;
basicTF.df = (basicTF.fmax-basicTF.fmin)/512 ;

if strcmp(basicTF.representation,'conceFT')
    basicTF.MT = 10 ;
end

VideoWriting = 0 ; % Change to 1 when recording

%% Real-time SST

[SSTtot, ForecastTime, TFRtime] = BoundEffRed_RT(xtot,fs,forecastMethod,basicTF,VideoWriting) ;

dtmax = basicTF.hop/fs ;
Hmin = ceil( fs/2 * ( mean(ForecastTime) + sqrt(mean(ForecastTime)^2 + 4*mean(TFRtime)/fs) ) );

%save('../../Results/RealTime','ForecastTime','TFRtime','dtmax')
