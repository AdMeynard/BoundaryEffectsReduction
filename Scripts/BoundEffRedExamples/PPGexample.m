clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load PPG signal
addpath ../../Signals/
data = load('PPG');

xtot = data.PPG;
fs = data.fs; % sampling frequency

x = xtot(30e3:34e3);
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

N = length(x);

%% Forecasting

HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets
if extK + extM >length(x) - 10 
    extK = extK/2 ; extM = extM/2 ;
end


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
tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*round(extSEC*fs)) ;

xxTRUE = ( xtot((30e3-L):(34e3+L)) - mu ) / s ;
xxZP = [ zeros(L,1); x; zeros(L,1) ] ; % zero padding


figure;
plot(tt,xxLSE,tt,xxEDMD,tt,xxGPR,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
set(gca,'fontsize',14);
legend('LSE Estimated Extended signal','EDMD Estimated Extended signal','GPR Estimated Extended signal','Ground truth Extended signal','Original signal'); 
xlabel('Time (s)'); ylabel('Signals'); title('Time series'); axis tight;

%% SST
basicTF.hop = 10;
basicTF.win = 1501;
fmin = 0 ;
fmax = 0.01;

% On the original signal
[Vx, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(double(x-mean(x)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[VxxTRUE, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the zero-padded extended signal 
[VxxZP, SSTxxZP, ~, tfrsqticZP] = ConceFT_sqSTFT_C(double(xxZP-mean(xxZP)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxLSE, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxEDMD, SSTxxEDMD, ~, ~] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxGPR, SSTxxGPR, ~, ~] = ConceFT_sqSTFT_C(double(xxGPR-mean(xxGPR)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

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

%save('results','tfrsq3','tfrsq3EXT');
%% Optimal transport distance SST
tmp = tt(1:basicTF.hop:end) ;
tmp = (tmp>=0) & (tmp<=t(end)) ;

tSSTx = t(1:basicTF.hop:end);
tSSText = tt(1:basicTF.hop:end);
n0 = find(tSSText>=min(tSSTx),1);
nf = n0 + length(tSSTx)-1;

tfrsqTRUEw = SSTxxTRUE(:,n0:nf) ;
tfrsqZPw = SSTxxZP(:,n0:nf) ;
tfrsqLSEw = SSTxxLSE(:,n0:nf) ;
tfrsqEDMDw = SSTxxEDMD(:,n0:nf) ;
tfrsqGPRw = SSTxxGPR(:,n0:nf) ;

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

%% Reassignment
basicTF.hop = 10;
basicTF.win = 1501;
fmin = 0 ;
fmax = 0.01;

% On the original signal
[~, RSx, ~, tfrticx] = ConceFT_rsSTFT(double(x-mean(x)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the true extended signal 
[~, RSxxTRUE, ~, tfrticxx] = ConceFT_rsSTFT(double(xxTRUE-mean(xxTRUE)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the zero-padded extended signal 
[~, RSxxZP] = ConceFT_rsSTFT(double(xxZP-mean(xxZP)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the estimated extended signal 
[~, RSxxLSE] = ConceFT_rsSTFT(double(xxLSE-mean(xxLSE)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the estimated extended signal 
[~, RSxxEDMD] = ConceFT_rsSTFT(double(xxEDMD-mean(xxEDMD)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
% On the estimated extended signal 
[~, RSxxGPR] = ConceFT_rsSTFT(double(xxGPR-mean(xxGPR)), fmin, fmax,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;

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

%% Optimal transport distance Reassignment
tmp = tt ;
tmp = (tmp>=0) & (tmp<=t(end)) ;

rss = abs(RSx).^2 ;
rssTRUEw = abs(RSxxTRUE(:,tmp)).^2 ;
rssZPw = abs(RSxxZP(:,tmp)).^2 ;
rssLSEw = abs(RSxxLSE(:,tmp)).^2 ;
rssEDMDw = abs(RSxxEDMD(:,tmp)).^2 ;
rssGPRw = abs(RSxxGPR(:,tmp)).^2 ;

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


%% STFT
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

%% Optimal transport distance STFT

stft = abs(Vx).^2 ;
stftTRUEw = abs(VxxTRUE(:,n0:nf)).^2 ;
stftZPw = abs(VxxZP(:,n0:nf)).^2 ;
stftLSEw = abs(VxxLSE(:,n0:nf)).^2 ;
stftEDMDw = abs(VxxEDMD(:,n0:nf)).^2 ;
stftGPRw = abs(VxxGPR(:,n0:nf)).^2 ;

stftOTDshort = slicedOT(stft, stftTRUEw) ;
stftOTDZP = slicedOT(stftZPw, stftTRUEw) ;
stftOTDLSE = slicedOT(stftLSEw, stftTRUEw) ;
stftOTDEDMD = slicedOT(stftEDMDw, stftTRUEw) ;
stftOTDGPR = slicedOT(stftGPRw, stftTRUEw) ;


fprintf('==============STFT===============\n')
fprintf(' _____________________________\n')
fprintf('|  Extension Method | Index D | \n')
fprintf('|-------------------|---------|\n')
fprintf('|   Zero-padding    |  %.3f  |\n', stftOTDZP/stftOTDshort)
fprintf('|      SigExt       |  %.3f  |\n', stftOTDLSE/stftOTDshort)
fprintf('|   EDMD extension  |  %.3f  |\n', stftOTDEDMD/stftOTDshort)
fprintf('|   GPR extension   |  %.3f  |\n', stftOTDGPR/stftOTDshort)
fprintf('|___________________|_________|\n')