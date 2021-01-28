clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load Arterial Blood Pressure signal
addpath ../../Signals/
data = load('ABP');

xtot = data.ABP;
fs = data.fs ;  % sampling frequency

HOP0 = 8 ;
xtot = xtot(1:HOP0:end) ;
fs = fs/HOP0 ;

n0 = 1 ;
nf = 100*fs+1 ;
N = nf-n0+1 ;

x = xtot(n0:nf);
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

xf = [xtot((nf+1):end); xtot(1:nf)] ;
xf = (xf - mu)/s ;

%% Forecasting

HOP = 1 ;
extSEC = 9.2 ; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets

method.name = 'SigExt' ;
tic; 
xxLSE = forecasting(xf,L,HOP,extK,extM,method);
xxLSE = [x; xxLSE] ;
LSEtime = toc;

xxTRUE = ( xtot(n0:(nf+L)) - mu ) / s ;

t = linspace(0, (N-1)/fs, N) ;
tt = linspace(0, (N-1+L)/fs, N+L) ;

%% SST and STFT
basicTF.hop = 10 ;
alpha = 0.1*sqrt(2*log(100)) ; % ratio window width/window length
twin = 2*extSEC/alpha ;
basicTF.win = 2*floor(twin*fs/2)+1 ;
fmin = 0/fs ;
fmax = 1.5/fs ;
nfreq = 512 ;
df = (fmax-fmin)/nfreq ;

% On the original signal
[Vx, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(x, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[VxxTRUE, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(xxTRUE, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxLSE, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

tc = t-40 ;
ttc = tt-40 ;

figure;
subplot('Position',[0.055 0.12 0.435 0.865]);
imagesc(tc(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(SSTx)/5e1)); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
xlim([0 tc(end)]) ;
subplot('Position',[0.56 0.12 0.435 0.865]);
imagesc(ttc(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
xlim([0 tc(end)]) ;

dim = [.4625 0.12 .027 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;
dim = [.968 0.12 .027 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;

ax = [0.24 0.24] ;
ay = [0.28 0.22] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)
ax = [0.27 0.27] ;
ay = [0.15 0.21] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)
ax = [0.48 0.48] ;
ay = [0.15 0.21] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)
ax = [0.98 0.98] ;
ay = [0.14 0.20] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)

ax = [0.24 0.24] ;
ay = [0.705 0.645] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.27 0.27] ;
ay = [0.59 0.65] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.48 0.48] ;
ay = [0.59 0.65] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.98 0.98] ;
ay = [0.585 0.645] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)

figure;
subplot('Position',[0.055 0.12 0.435 0.865]);
imagesc(tc(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(Vx)/5e1)); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
xlim([0 tc(end)]) ;
subplot('Position',[0.56 0.12 0.435 0.865]);
imagesc(ttc(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(VxxLSE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
xlim([0 tc(end)]) ;

dim = [.4625 0.12 .027 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;
dim = [.968 0.12 .027 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;

%% Performance index for SST and STFT

tSSTx = t(1:basicTF.hop:end);
tSSText = tt(1:basicTF.hop:end);
n0 = find(tSSText>=min(tSSTx),1);
nf = n0 + length(tSSTx)-1;

stft = abs(Vx).^2 ;
stftTRUEw = abs(VxxTRUE(:,n0:nf)).^2 ;
stftLSEw = abs(VxxLSE(:,n0:nf)).^2 ;

tfrsqTRUEw = SSTxxTRUE(:,n0:nf) ;
tfrsqLSEw = SSTxxLSE(:,n0:nf) ;

stftOTDshort = slicedOT(stft, stftTRUEw) ;
stftOTDLSE = slicedOT(stftLSEw, stftTRUEw) ;

sstOTDshort = slicedOT(SSTx, tfrsqTRUEw) ;
sstOTDLSE = slicedOT(tfrsqLSEw, tfrsqTRUEw) ;

fprintf('=======TF Representations====\n')
fprintf(' ____________________________\n')
fprintf('| Extension Method | Index D | \n')
fprintf('|------------------|---------|\n')
fprintf('|       STFT       |  %6.3f |\n', stftOTDLSE/stftOTDshort)
fprintf('|        SST       |  %6.3f |\n', sstOTDLSE/sstOTDshort)
fprintf('|__________________|_________|\n')