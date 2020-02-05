clear all; close all; clc;
addpath(genpath('../SST'));
addpath('../algorithm/');

%% load THO signal
addpath ../Signals
data = load('THO.mat');
xtot = data.THO;

xtot = xtot(6.8e4:36.8e4) ; % avoid pulses and other weird patterns
fs = 100 ; % sampling frequency

N = 6e3; % subsignals length

% Forecasting parameters
HOP = 1 ;
extSEC = 7 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round(1.5*L) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

[forecastErrLSE, sstS, sstEXT, sstLSE] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC);

% figure;
% plot(forecastErrLSEV);

save('../Results/resultSucForTHO','forecastErrLSE','sstS','sstEXT','sstLSE');
%% Plot
% 
% %% SST
% basicTF.hop = 10;
% basicTF.win = 1501;
% 
% % On the original signal
% [~, ~, tfrsq3, ConceFT3, tfrsqtic] = ConceFT_sqSTFT_C(double(x-mean(x)), 0, 0.01,...
%             1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% % On the extended signal 
% [~, ~, tfrsq3EXT, ConceFT3EXT, tfrsqticEXT] = ConceFT_sqSTFT_C(double(xx-mean(xx)), 0, 0.01,...
%             1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% 
% figure;
% subplot(1,2,1);
% imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(tfrsq3)/8e5));
% xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original signal');
% subplot(1,2,2);
% imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(tfrsq3EXT)/8e5)); xlim([0 t(end)]);
% xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on extended signal');
