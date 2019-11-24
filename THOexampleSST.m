clear all; close all; clc;
addpath(genpath('SST'));

%% load THO signal
addpath NonStationaritySleepDataSet/1
data = load('THO.mat');
xtot = data.THO;
fs = 100 ; % sampling frequency

x = xtot(204e3:208e3);
N = length(x);

%% Forecasting

HOP = 1 ;
extSEC = 2 ; % the extension is of extSEC second
extP = 5*extSEC*fs ; % dimension of embedding / signals length
extN = 2*extP ;  % number of points to estimate A / size of datasets
if extN + extP >length(x) - 10 
    extN = extN/2 ; extP = extP/2 ;
end


xx = SigExtension(x,fs,HOP,extN,extP,extSEC);

%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;

xxTRUE = xtot((204e3-extSEC*fs):(208e3+extSEC*fs));

figure;
plot(tt,xx,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
legend('Estimated Extended signal','Ground truth Extended signal','Original signal'); 
xlabel('Time (s)'); ylabel('Signals'); title('Time series');

%% SST
basicTF.hop = 10;
basicTF.win = 1501;

% On the original signal
[~, ~, tfrsq3, ConceFT3, tfrsqtic] = ConceFT_sqSTFT_C(double(x-mean(x)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the extended signal 
[~, ~, tfrsq3EXT, ConceFT3EXT, tfrsqticEXT] = ConceFT_sqSTFT_C(double(xx-mean(xx)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

figure;
subplot(1,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(tfrsq3)/8e5));
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original signal');
subplot(1,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(tfrsq3EXT)/8e5)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on extended signal');

