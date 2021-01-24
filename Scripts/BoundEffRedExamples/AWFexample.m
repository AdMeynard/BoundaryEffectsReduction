clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load AWF
addpath ../../Signals/
data = load('AWF');

xtot = data.AWF ;
fs = data.fs ; % sampling frequency

HOP0 = 5 ;
xtot = xtot(1:HOP0:end) ;
fs = fs/HOP0 ;

n0 = 1 ;
N = ceil(4*length(xtot)/5) ; %10e3 ;
nf = n0 + N-1 ;

x = xtot(n0:nf);
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

%% Forecasting

HOP = 1 ;
extSEC = 57 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 1.8*extM );  % number of points to estimate A / size of datasets

method.name = 'SigExt' ;
tic; 
xxLSE = forecasting(x,L,HOP,extK,extM,method);
xxLSE = [x; xxLSE] ;
LSEtime = toc;

% method.name = 'edmd' ;
% method.param = 100 ;
% tic; 
% xxEDMD = forecasting(x,L,HOP,extK,extM,method); 
% xxEDMD = [x; xxEDMD] ;
% EDMDtime = toc;
% 
% method.name = 'gpr' ;
% tic; 
% xxGPR = forecasting(x,L,HOP,extK,extM,method); 
% xxGPR = [x; xxGPR] ;
% GPRtime = toc;
% 
% fprintf(' __________________________________________\n')
% fprintf('| Extension Method | Computing time (sec.) |\n')
% fprintf('|------------------------------------------|\n')
% fprintf('|     SigExt       |         %.2f          |\n', LSEtime)
% fprintf('|  EDMD extension  |         %.2f          |\n', EDMDtime)
% fprintf('|  GPR extension   |       %.2f          |\n', GPRtime)
% fprintf('|__________________|_______________________|\n')

%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(0, (N-1+L)/fs, N+L) ;

xxTRUE = ( xtot(n0:(nf+L)) - mu ) / s ;
xxZP = [ zeros(L,1); x; zeros(L,1) ] ; % zero padding

figure;

% subplot(3,1,1)
plot(tt,xxLSE,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
ylabel('Signals'); axis tight; xlim([t(end)-1.5*extSEC tt(end)]);
xticks([]);
legend({'{\sf SigExt} extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
set(gca,'fontsize',18);


% subplot(3,1,2)
% plot(tt,xxEDMD,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
% ylabel('Signals'); axis tight; xlim([t(end)-1.5*extSEC tt(end)]);
% xticks([]);
% legend({'EDMD extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
% set(gca,'fontsize',18);
% 
% subplot(3,1,3)
% plot(tt,xxGPR,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
% xlabel('Time (s)'); ylabel('Signals'); axis tight; xlim([t(end)-1.5*extSEC tt(end)]);
% legend({'GPR extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
% set(gca,'fontsize',18);

%% SST and STFT
basicTF.hop = 5 ;
twin = 114 ;
basicTF.win = 2*floor(twin*fs/2)+1 ;
fmin = 0/fs ;
fmax = 2/fs ;
nfreq = 1024 ;
df = (fmax-fmin)/nfreq ;

% On the original signal
[Vx, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(x, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[VxxTRUE, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(xxTRUE, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the zero-padded extended signal 
[VxxZP, SSTxxZP, ~, tfrsqticZP] = ConceFT_sqSTFT_C(xxZP, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxLSE, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% % On the estimated extended signal 
% [VxxEDMD, SSTxxEDMD, ~, ~] = ConceFT_sqSTFT_C(xxEDMD, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% % On the estimated extended signal 
% [VxxGPR, SSTxxGPR, ~, ~] = ConceFT_sqSTFT_C(xxGPR, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(SSTx)/5e1)); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(SSTxxTRUE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on LSE estimated long signal (short signal extended)');
% subplot(2,2,4);
% imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxEDMD)/1e2)); xlim([0 t(end)]); axis xy ;
% xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on EDMD estimated long signal (short signal extended)');

% figure;
% imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
% xlabel('Time (s)'); ylabel('Frequency (Hz)');
% set(gca,'fontsize',22) ;

figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(Vx)/1e2)); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(VxxTRUE)/1e2)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(VxxLSE)/1e2)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT LSE estimated long signal (short signal extended)');
% subplot(2,2,4);
% imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(VxxEDMD)/1e2)); xlim([0 t(end)]); axis xy ;
% xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('STFT on EDMD estimated long signal (short signal extended)');

%save('results','tfrsq3','tfrsq3EXT');
%% Performance index for SST and STFT

tSSTx = t(1:basicTF.hop:end);
tSSText = tt(1:basicTF.hop:end);
n0 = find(tSSText>=min(tSSTx),1,'first') ;
nf = find(tSSText<=max(tSSTx),1,'last') ;

stft = abs(Vx).^2 ;
stftTRUEw = abs(VxxTRUE(:,n0:nf)).^2 ;
stftZPw = abs(VxxZP(:,n0:nf)).^2 ;
stftLSEw = abs(VxxLSE(:,n0:nf)).^2 ;
% stftEDMDw = abs(VxxEDMD(:,n0:nf)).^2 ;
% stftGPRw = abs(VxxGPR(:,n0:nf)).^2 ;

tfrsqTRUEw = SSTxxTRUE(:,n0:nf) ;
tfrsqZPw = SSTxxZP(:,n0:nf) ;
tfrsqLSEw = SSTxxLSE(:,n0:nf) ;
% tfrsqEDMDw = SSTxxEDMD(:,n0:nf) ;
% tfrsqGPRw = SSTxxGPR(:,n0:nf) ;

stftOTDshort = slicedOT(stft, stftTRUEw) ;
stftOTDZP = slicedOT(stftZPw, stftTRUEw) ;
stftOTDLSE = slicedOT(stftLSEw, stftTRUEw) ;
% stftOTDEDMD = slicedOT(stftEDMDw, stftTRUEw) ;
% stftOTDGPR = slicedOT(stftGPRw, stftTRUEw) ;

sstOTDshort = slicedOT(SSTx, tfrsqTRUEw) ;
sstOTDZP = slicedOT(tfrsqZPw, tfrsqTRUEw) ;
sstOTDLSE = slicedOT(tfrsqLSEw, tfrsqTRUEw) ;
% sstOTDEDMD = slicedOT(tfrsqEDMDw, tfrsqTRUEw) ;
% sstOTDGPR = slicedOT(tfrsqGPRw, tfrsqTRUEw) ;


fprintf('=============SST==============\n')
fprintf(' _____________________________\n')
fprintf('|  Extension Method | Index D | \n')
fprintf('|-------------------|---------|\n')
fprintf('|   Zero-padding    |  %6.3f |\n', sstOTDZP/sstOTDshort)
fprintf('|      SigExt       |  %6.3f |\n', sstOTDLSE/sstOTDshort)
% fprintf('|   EDMD extension  |  %6.3f |\n', sstOTDEDMD/sstOTDshort)
% fprintf('|   GPR extension   |  %6.3f |\n', sstOTDGPR/sstOTDshort)
fprintf('|___________________|_________|\n')

fprintf('=============STFT==============\n')
fprintf(' _____________________________\n')
fprintf('|  Extension Method | Index D | \n')
fprintf('|-------------------|---------|\n')
fprintf('|   Zero-padding    |  %6.3f |\n', stftOTDZP/stftOTDshort)
fprintf('|      SigExt       |  %6.3f |\n', stftOTDLSE/stftOTDshort)
% fprintf('|   EDMD extension  |  %5.3f  |\n', stftOTDEDMD/stftOTDshort)
% fprintf('|   GPR extension   |  %5.3f  |\n', stftOTDGPR/stftOTDshort)
fprintf('|___________________|_________|\n')