clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load PPG signal
addpath ../../Signals/
data = load('ECG');

xtot = data.ECG ;
fs = 200 ; % sampling frequency

n0 = 51e3;
nf = 58e3;
x = xtot(n0:nf);
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

N = length(x);

%% Forecasting

HOP = 1 ;
extSEC = 6 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets

method.name = 'SigExt' ;
tic; 
xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
LSEtime = toc;

method.name = 'edmd' ;
method.param = 100 ;
tic; 
xxEDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
EDMDtime = toc;

method.name = 'gpr' ;
tic; 
xxGPR = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
GPRtime = toc;

fprintf(' __________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) |\n')
fprintf('|------------------|-----------------------|\n')
fprintf('|     SigExt       |        %.2f          |\n', LSEtime)
fprintf('|  EDMD extension  |        %.2f          |\n', EDMDtime)
fprintf('|  GPR extension   |       %.2f          |\n', GPRtime)
fprintf('|__________________|_______________________|\n')


%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(-L/fs, (N-1+L)/fs, N+2*L) ;

xxTRUE = ( xtot((n0-L):(nf+L)) - mu ) / s ;
xxZP = [ zeros(L,1); x; zeros(L,1) ] ; % zero padding


figure;

subplot(3,1,1)
plot(tt,xxLSE,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
ylabel('Signals'); axis tight; xlim([t(end)-1.5*extSEC tt(end)]);
xticks([]);
legend({'{\sf SigExt} extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
set(gca,'fontsize',18);

subplot(3,1,2)
plot(tt,xxEDMD,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
ylabel('Signals'); axis tight; xlim([t(end)-1.5*extSEC tt(end)]);
xticks([]);
legend({'EDMD extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
set(gca,'fontsize',18);

subplot(3,1,3)
plot(tt,xxGPR,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
xlabel('Time (s)'); ylabel('Signals'); axis tight; xlim([t(end)-1.5*extSEC tt(end)]);
legend({'GPR extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
set(gca,'fontsize',18);

%% SST
basicTF.hop = 10 ;
basicTF.win = 2*L+1 ;
basicTF.fmin = 0 ;
basicTF.fmax = 0.02 ;
basicTF.df = 1e-5 ;

% On the original signal
[~, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(x-mean(x), basicTF.fmin, basicTF.fmax,...
            basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[~, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(xxTRUE, basicTF.fmin, basicTF.fmax,...
            basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the zero-padded extended signal 
[~, SSTxxZP, ~, tfrsqticZP] = ConceFT_sqSTFT_C(xxZP, basicTF.fmin, basicTF.fmax,...
            basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE, basicTF.fmin, basicTF.fmax,...
            basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, SSTxxEDMD, ~, ~] = ConceFT_sqSTFT_C(xxEDMD, basicTF.fmin, basicTF.fmax,...
            basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, SSTxxGPR, ~, ~] = ConceFT_sqSTFT_C(xxGPR, basicTF.fmin, basicTF.fmax,...
            basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(SSTx)/1e2));
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(SSTxxTRUE)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on estimated long signal (short signal extended)');
subplot(2,2,4);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxEDMD)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on estimated long signal (short signal extended)');

%save('results','tfrsq3','tfrsq3EXT');
%% Optimal transport distance
tmp = tt(1:basicTF.hop:end) ;
tmp = (tmp>=0) & (tmp<=t(end)) ;

tfrsqTRUEw = SSTxxTRUE(:,tmp) ;
tfrsqZPw = SSTxxZP(:,tmp) ;
tfrsqLSEw = SSTxxLSE(:,tmp) ;
tfrsqEDMDw = SSTxxEDMD(:,tmp) ;
tfrsqGPRw = SSTxxGPR(:,tmp) ;

OTDshort = slicedOT(SSTx, tfrsqTRUEw) ;
OTDZP = slicedOT(tfrsqZPw, tfrsqTRUEw) ;
OTDLSE = slicedOT(tfrsqLSEw, tfrsqTRUEw) ;
OTDEDMD = slicedOT(tfrsqEDMDw, tfrsqTRUEw) ;
OTDGPR = slicedOT(tfrsqGPRw, tfrsqTRUEw) ;

fprintf('=============SST===============\n')
fprintf(' ____________________________\n')
fprintf('| Extension Method | Index D | \n')
fprintf('|------------------|---------|\n')
fprintf('|     SigExt       |  %.3f  |\n', OTDLSE/OTDshort)
fprintf('|  EDMD extension  |  %.3f  |\n', OTDEDMD/OTDshort)
fprintf('|  GPR extension   |  %.3f  |\n', OTDGPR/OTDshort)
fprintf('|__________________|_________|\n')
