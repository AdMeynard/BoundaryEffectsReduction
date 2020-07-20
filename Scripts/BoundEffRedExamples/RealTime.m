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
fs = 100/subT ; % sampling frequency

minutes = 3 ; 
Nmax = minutes*60*fs ;
xtot = xtot(1:Nmax) ; % shorten signal

%% Parameters
forecastMethod.name = 'lse' ;

basicTF.representation = 'conceFT' ;

basicTF.hop = round(20/subT) ;
basicTF.win = 501 ; % window length (in samples)
basicTF.fmin = 0.5/fs ;
basicTF.fmax = 4/fs ;
basicTF.df = 2e-4 ;

if strcmp(basicTF.representation,'conceFT')
    basicTF.MT = 10 ;
end

VideoWriting = 0 ; % Change to 1 when recording

%% Real-time SST

[SSTtot, dt] = BoundEffRed_RT(xtot,fs,forecastMethod,basicTF,VideoWriting) ;

dtmax = basicTF.hop/fs ;
save('../../Results/RealTime','dt','dtmax')
