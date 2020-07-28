%% Fig PPG
clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load PPG signal
addpath ../../Signals
data = load('PPG.mat');
xtot = data.PPG;

fs = data.fs ; % sampling frequency

n0 = 6e3 ;
N = 4e3; % subsignals length
x0 = xtot((n0):(n0+N-1)) ; % Subsignal
mu = mean(x0);
sigma = std(x0);
x0 = (x0 - mu) / sigma ;

dgamma = cumsum(randn(N,1)) ;
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
dgamma = filter(b,a,dgamma);
dgamma = dgamma - mean(dgamma) ;
dgamma = 0.2*dgamma/max(abs(dgamma)) ;
dgamma = dgamma + 1 ;

gamma = cumsum(dgamma) - dgamma(1);
gamma = gamma/gamma(end) ;
x = sigwarp( x0,gamma,dgamma ) ;

%% Extension
HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets
method.name = 'SigExt' ;

xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ;

%% SST
basicTF.hop = 10 ;
basicTF.win = 2*L+1 ;
basicTF.fmin = 0/fs ;
basicTF.fmax = 5.5/fs ;
basicTF.df = 1e-5 ;

t = linspace(0, (N-1)/fs, N) ;
tsst = t(1:basicTF.hop:end);

[~, sstS, ~,fsst] = ConceFT_sqSTFT_C(x, basicTF.fmin, basicTF.fmax,basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
fsst = fs*fsst ;

tt = linspace(-L/fs, (N-1+L)/fs, N+2*L) ;
tsstEXT = tt(1:basicTF.hop:end);
[~, sstLSE, ~, ~] = ConceFT_sqSTFT_C(xxLSE-mean(xxLSE), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;


figure;
subplot(3,2,[1 2]);
plot(t,x,'k','linewidth',2); axis tight ; xlim([t(end)-15 t(end)]) ;
xlabel('Time (s)'); ylabel('Signal'); 
set(gca,'fontsize',18); grid on;

subplot(3,2,[3 5]);
imagesc(tsst,fsst,log1p(abs(sstS)/1e1)); axis xy; xlim([t(end)-15 t(end)])
xlabel('Time (s)'); ylabel('Frequency (Hz)'); colormap(1-gray)
set(gca,'fontsize',18);

subplot(3,2,[4 6]);
imagesc(tsstEXT,fsst,log1p(abs(sstLSE)/1e1)); axis xy; xlim([t(end)-15 t(end)])
xlabel('Time (s)'); ylabel('Frequency (Hz)'); colormap(1-gray)
set(gca,'fontsize',18);