clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load PPG signal
addpath ../../Signals/
data = load('PPG');

xtot = data.PPG;
fs = data.fs; % sampling frequency

n0 = 30e3 ;
N = 4e3 ;
nf = n0 + N-1 ;

x = xtot(n0:nf);
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

%% Forecasting

HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
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
fprintf('|------------------------------------------|\n')
fprintf('|     SigExt       |         %.2f          |\n', LSEtime)
fprintf('|  EDMD extension  |         %.2f          |\n', EDMDtime)
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

%% SST and STFT
basicTF.hop = 10;
basicTF.win = 1501;
fmin = 0 ;
fmax = 0.01 ;
df = 1e-5 ;

% On the original signal
[Vx, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(x, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[VxxTRUE, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(xxTRUE, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the zero-padded extended signal 
[VxxZP, SSTxxZP, ~, tfrsqticZP] = ConceFT_sqSTFT_C(xxZP, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxLSE, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxEDMD, SSTxxEDMD, ~, ~] = ConceFT_sqSTFT_C(xxEDMD, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxGPR, SSTxxGPR, ~, ~] = ConceFT_sqSTFT_C(xxGPR, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(SSTx)/1e2));
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(SSTxxTRUE)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on LSE estimated long signal (short signal extended)');
subplot(2,2,4);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxEDMD)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on EDMD estimated long signal (short signal extended)');

figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(Vx)/1e2)); ylim([0 10]) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(VxxTRUE)/1e2)); xlim([0 t(end)]); ylim([0 10]) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(VxxLSE)/1e2)); xlim([0 t(end)]); ylim([0 10]) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT LSE estimated long signal (short signal extended)');
subplot(2,2,4);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(VxxEDMD)/1e2)); xlim([0 t(end)]); ylim([0 10]) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT on EDMD estimated long signal (short signal extended)');

%save('results','tfrsq3','tfrsq3EXT');
%% Performance index for SST and STFT

tSSTx = t(1:basicTF.hop:end);
tSSText = tt(1:basicTF.hop:end);
n0 = find(tSSText>=min(tSSTx),1);
nf = n0 + length(tSSTx)-1;

stft = abs(Vx).^2 ;
stftTRUEw = abs(VxxTRUE(:,n0:nf)).^2 ;
stftZPw = abs(VxxZP(:,n0:nf)).^2 ;
stftLSEw = abs(VxxLSE(:,n0:nf)).^2 ;
stftEDMDw = abs(VxxEDMD(:,n0:nf)).^2 ;
stftGPRw = abs(VxxGPR(:,n0:nf)).^2 ;

tfrsqTRUEw = SSTxxTRUE(:,n0:nf) ;
tfrsqZPw = SSTxxZP(:,n0:nf) ;
tfrsqLSEw = SSTxxLSE(:,n0:nf) ;
tfrsqEDMDw = SSTxxEDMD(:,n0:nf) ;
tfrsqGPRw = SSTxxGPR(:,n0:nf) ;

stftOTDshort = slicedOT(stft, stftTRUEw) ;
stftOTDZP = slicedOT(stftZPw, stftTRUEw) ;
stftOTDLSE = slicedOT(stftLSEw, stftTRUEw) ;
stftOTDEDMD = slicedOT(stftEDMDw, stftTRUEw) ;
stftOTDGPR = slicedOT(stftGPRw, stftTRUEw) ;

sstOTDshort = slicedOT(SSTx, tfrsqTRUEw) ;
sstOTDZP = slicedOT(tfrsqZPw, tfrsqTRUEw) ;
sstOTDLSE = slicedOT(tfrsqLSEw, tfrsqTRUEw) ;
sstOTDEDMD = slicedOT(tfrsqEDMDw, tfrsqTRUEw) ;
sstOTDGPR = slicedOT(tfrsqGPRw, tfrsqTRUEw) ;


fprintf('==============SST===============\n')
fprintf(' _______________________________\n')
fprintf('|  Extension Method |  Index D  | \n')
fprintf('|-------------------|-----------|\n')
fprintf('|   Zero-padding    | %.3f |\n', sstOTDZP/sstOTDshort)
fprintf('|      SigExt       | %.3f |\n', sstOTDLSE/sstOTDshort)
fprintf('|   EDMD extension  | %.3f |\n', sstOTDEDMD/sstOTDshort)
fprintf('|   GPR extension   | %.3f |\n', sstOTDGPR/sstOTDshort)
fprintf('|___________________|___________|\n')

fprintf('==============STFT===============\n')
fprintf(' _____________________________\n')
fprintf('|  Extension Method | Index D | \n')
fprintf('|-------------------|---------|\n')
fprintf('|   Zero-padding    |  %.3f  |\n', stftOTDZP/stftOTDshort)
fprintf('|      SigExt       |  %.3f  |\n', stftOTDLSE/stftOTDshort)
fprintf('|   EDMD extension  |  %.3f  |\n', stftOTDEDMD/stftOTDshort)
fprintf('|   GPR extension   |  %.3f  |\n', stftOTDGPR/stftOTDshort)
fprintf('|___________________|_________|\n')

%% Reassignment

% On the original signal
[~, RSx, ~, tfrticx] = ConceFT_rsSTFT(x, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the true extended signal 
[~, RSxxTRUE, ~, tfrticxx] = ConceFT_rsSTFT(xxTRUE, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the zero-padded extended signal 
[~, RSxxZP] = ConceFT_rsSTFT(xxZP, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the estimated extended signal 
[~, RSxxLSE] = ConceFT_rsSTFT(xxLSE, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the estimated extended signal 
[~, RSxxEDMD] = ConceFT_rsSTFT(xxEDMD, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the estimated extended signal 
[~, RSxxGPR] = ConceFT_rsSTFT(xxGPR, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1) ;

figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrticx*fs,log1p(abs(RSx)/1e2));
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('RS on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrticxx*fs,log1p(abs(RSxxTRUE)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('RS on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrticxx*fs,log1p(abs(RSxxLSE)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('RS on LSE estimated long signal (short signal extended)');
subplot(2,2,4);
imagesc(tt(1:basicTF.hop:end),tfrticxx*fs,log1p(abs(RSxxEDMD)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('RS on EDMD estimated long signal (short signal extended)');

%save('results','tfrsq3','tfrsq3EXT');

%% Performance index for Reassignment

rss = abs(RSx).^2 ;
rssTRUEw = abs(RSxxTRUE(:,n0:nf)).^2 ;
rssZPw = abs(RSxxZP(:,n0:nf)).^2 ;
rssLSEw = abs(RSxxLSE(:,n0:nf)).^2 ;
rssEDMDw = abs(RSxxEDMD(:,n0:nf)).^2 ;
rssGPRw = abs(RSxxGPR(:,n0:nf)).^2 ;

rssOTDshort = slicedOT(rss, rssTRUEw) ;
rssOTDZP = slicedOT(rssZPw, rssTRUEw) ;
rssOTDLSE = slicedOT(rssLSEw, rssTRUEw) ;
rssOTDEDMD = slicedOT(rssEDMDw, rssTRUEw) ;
rssOTDGPR = slicedOT(rssGPRw, rssTRUEw) ;

fprintf('=================RS==============\n')
fprintf(' ________________________________\n')
fprintf('|  Extension Method  |  Index D  | \n')
fprintf('|--------------------|-----------|\n')
fprintf('|      SigExt        | %.3f |\n', rssOTDLSE/rssOTDshort)
fprintf('|  EDMD extension    | %.3f |\n', rssOTDEDMD/rssOTDshort)
fprintf('|  GPR extension     | %.3f |\n', rssOTDGPR/rssOTDshort)
fprintf('|____________________|___________|\n')

%% CWT
wav_typ = 'sharp';
Ms = 300 ;
wav_par = 9 ;
fmax = 10 ;
fmin = fmax/10 ; 

figure;
% On the original signal
subplot(2,2,1);
[WxxZP,~] = display_cwt_JEFAS(xxZP,fs,fmin,fmax,Ms,wav_typ,wav_par) ;
[Wx,~] = display_cwt_JEFAS(x,fs,fmin,fmax,Ms,wav_typ,wav_par) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('CWT on original short signal');
% On the true extended signal
subplot(2,2,2);
[WxxTRUE,~] = display_cwt_JEFAS(xxTRUE,fs,fmin,fmax,Ms,wav_typ,wav_par) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('CWT on original long signal');
% On the estimated extended signal
subplot(2,2,3);
[WxxLSE,~] = display_cwt_JEFAS(xxLSE,fs,fmin,fmax,Ms,wav_typ,wav_par) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('CWT on estimated long signal (short signal extended)');
% On the estimated extended signal
subplot(2,2,4);
[WxxGPR,~] = display_cwt_JEFAS(xxGPR,fs,fmin,fmax,Ms,wav_typ,wav_par) ;
[WxxEDMD,~] = display_cwt_JEFAS(xxEDMD,fs,fmin,fmax,Ms,wav_typ,wav_par) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('CWT on estimated long signal (short signal extended)');

%save('results','tfrsq3','tfrsq3EXT');

%% Optimal transport distance CWT
tmp = tt ;
tmp = (tmp>=0) & (tmp<=t(end)) ;

scalo = abs(Wx).^2 ;
scaloTRUEw = abs(WxxTRUE(:,tmp)).^2 ;
scaloZPw = abs(WxxZP(:,tmp)).^2 ;
scaloLSEw = abs(WxxLSE(:,tmp)).^2 ;
scaloEDMDw = abs(WxxEDMD(:,tmp)).^2 ;
scaloGPRw = abs(WxxGPR(:,tmp)).^2 ;

cwtOTDshort = slicedOT(scalo, scaloTRUEw) ;
cwtOTDZP = slicedOT(scaloZPw, scaloTRUEw) ;
cwtOTDLSE = slicedOT(scaloLSEw, scaloTRUEw) ;
cwtOTDEDMD = slicedOT(scaloEDMDw, scaloTRUEw) ;
cwtOTDGPR = slicedOT(scaloGPRw, scaloTRUEw) ;


fprintf('================CWT==============\n')
fprintf(' ________________________________\n')
fprintf('|  Extension Method  |  Index D  | \n')
fprintf('|--------------------|-----------|\n')
fprintf('|   Zero-padding     | %.3f |\n', cwtOTDZP/cwtOTDshort)
fprintf('|      SigExt        | %.3f |\n', cwtOTDLSE/cwtOTDshort)
fprintf('|   EDMD extension   | %.3f |\n', cwtOTDEDMD/cwtOTDshort)
fprintf('|   GPR extension    | %.3f |\n', cwtOTDGPR/cwtOTDshort)
fprintf('|____________________|___________|\n')