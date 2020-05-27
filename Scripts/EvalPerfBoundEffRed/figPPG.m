%% Fig PPG
clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load ECG signal
addpath ../../Signals
data = load('PPG.mat');
xtot = data.PPG;

fs = data.fs ; % sampling frequency
N = 4e3; % subsignals length
% Subsignal
k = 5 ;
x = xtot((k*N+1):((k+1)*N)) ;
mu = mean(x);
sigma = std(x);
x = (x - mu) / sigma ;

%% SST
basicTF.hop = 10 ;
basicTF.win = 1501 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.06 ;

t = linspace(0, (N-1)/fs, N) ;
tsst = t(1:basicTF.hop:end);

[~, ~, sstS, ~] = ConceFT_sqSTFT_C(x, basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
ff = linspace(0,fs*basicTF.fmax,400) ;

figure;
imagesc(tsst,ff,log1p(abs(sstS)/1e1)); axis xy; 
xlabel('Time (sec.)'); ylabel('Frequency (Hz)'); colormap(1-gray);
set(gca,'fontsize',18);

figure;
imagesc(tsst,ff,log1p(abs(sstS)/1e1)); axis xy; xlim([27 t(end)])
xlabel('Time (sec.)'); ylabel('Frequency (Hz)'); colormap(1-gray)
set(gca,'fontsize',18);

%% Extension
HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets
method.name = 'lse' ;

xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ;

tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
tsstEXT = tt(1:basicTF.hop:end);
[~, ~, sstLSE, ~, ~] = ConceFT_sqSTFT_C(xxLSE-mean(xxLSE), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

figure;
imagesc(tsstEXT,ff,log1p(abs(sstLSE)/1e1)); axis xy; xlim([0 t(end)])
xlabel('Time (sec.)'); ylabel('Frequency (Hz)'); colormap(1-gray);
set(gca,'fontsize',18);

figure;
imagesc(tsstEXT,ff,log1p(abs(sstLSE)/1e1)); axis xy; xlim([27 t(end)])
xlabel('Time (sec.)'); ylabel('Frequency (Hz)'); colormap(1-gray)
set(gca,'fontsize',18);