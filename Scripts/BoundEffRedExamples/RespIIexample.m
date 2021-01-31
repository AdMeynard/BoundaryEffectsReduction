%% Run the extension methods, and the associated BoundEffRed TF representations on a 2-component respiratory signal
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load Respiratory signal
addpath ../../Signals/
data = load('RespII');

xtot = data.RESP ;
fs = 1/(data.timeRESP(2)-data.timeRESP(1)) ; % sampling frequency

HOP0 = 1 ;
xtot = xtot(1:HOP0:end) ;
fs = fs/HOP0 ;
tRESP = data.timeRESP(1:HOP0:end) ;

n0 = 1 ;
N = ceil(4*length(xtot)/5) ;
nf = n0 + N-1 ;

x = xtot(n0:nf) ;
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

%% Forecasting

HOP = 1 ;
extSEC = 10; % the extension is of extSEC second
L = round( extSEC*fs ) ;
extM = round( 1.5*L ) ; % dimension of embedding / signals length
extK = round( 2.5*extM );  % number of points to estimate A / size of datasets

method.name = 'SigExt' ;
tic; 
xxLSE = forecasting(x,L,HOP,extK,extM,method);
xxLSE = [x; xxLSE] ;
LSEtime = toc;

xxTRUE = ( xtot(n0:(nf+L)) - mu ) / s ;

t = tRESP(n0:nf) ; %linspace(0, (N-1)/fs, N) ;
tt = tRESP(n0:(nf+L)) ;%linspace(0, (N-1+L)/fs, N+L) ;

% display the SigExt extension
tm = t(end)-12.5 ;

figure;
subplot('Position',[0.045 0.52 0.942 0.465]);
plot(tt-tm,xxLSE,tt-tm,xxTRUE,'--',t-tm,x,'linewidth',2); grid on;
ylabel('Respiratory signal'); axis tight; xlim([0 t(end)-tm+7.5]);
xticklabels([]);
legend({'SigExt extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
set(gca,'fontsize',24);

xECG = data.ECG ;
xECG = (xECG - mean(xECG))/std(xECG) ;
tECG = data.timeECG ;

subplot('Position',[0.045 0.255 0.942 0.265]);
plot(tECG-tm,xECG,'k','linewidth',2); grid on;
xlabel('Time (s)'); ylabel('ECG signal'); axis tight; xlim([0 t(end)-tm+7.5]);
set(gca,'fontsize',24);

annotation('line',[0.142 0.142],[0.672 0.505],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.173 0.173],[0.584 0.502],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.205 0.205],[0.582 0.497],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.301 0.301],[0.661 0.506],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.333 0.333],[0.595 0.504],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.456 0.456],[0.595 0.477],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.487 0.487],[0.595 0.513],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.586 0.586],[0.623 0.504],'Color','m','LineWidth',3,'LineStyle',':');
annotation('line',[0.619 0.619],[0.603 0.516],'Color','m','LineWidth',3,'LineStyle',':');

%% SST and STFT
basicTF.hop = 3 ;
alpha = 0.1*sqrt(2*log(100)) ; % ratio window width/window length
twin = 2*extSEC/alpha ;
basicTF.win = 2*floor(twin*fs/2)+1 ;
fmin = 0/fs ;
fmax = 2/fs ;
nfreq = 1024 ;
df = (fmax-fmin)/nfreq ;

% On the original signal
[Vx, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(x, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[VxxTRUE, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(xxTRUE, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[VxxLSE, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

tc = t-10 ;
ttc = tt-10 ;

% SST
figure;
subplot('Position',[0.055 0.12 0.435 0.865]);
imagesc(tc(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(SSTx)/5e1)); axis xy ; colormap(1-gray);
xlim([0 tc(end)]) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
subplot('Position',[0.56 0.12 0.435 0.865]);
imagesc(ttc(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlim([0 tc(end)]) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;

dim = [.4325 0.12 .0575 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;
dim = [.938 0.12 .0575 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;

ax = [0.255 0.255] ;
ay = [0.325 0.265] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)
ax = [0.322 0.322] ;
ay = [0.235 0.295] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)
ax = [0.476 0.476] ;
ay = [0.212 0.272] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)
ax = [0.978 0.978] ;
ay = [0.212 0.272] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)

ax = [0.246 0.246] ;
ay = [0.830 0.770] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.344 0.344] ;
ay = [0.708 0.768] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.481 0.481] ;
ay = [0.709 0.769] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.971 0.971] ;
ay = [0.709 0.769] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)


figure;
subplot('Position',[0.062 0.516 0.9 0.37]);
imagesc(tfrsqtic*fs,tc(1:basicTF.hop:end),log1p(abs(SSTx')/5e1)); colormap(1-gray);
ylim([tc(end)-0.6*twin/2 tc(end)]); %yticks([200 220 240]);
ylabel('Time (s)','rotation',270,'VerticalAlignment','top'); xlabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
set(gca,'yticklabelrotation',270) ;
set(gca,'xaxisLocation','top')
subplot('Position',[0.062 0.02 0.9 0.37]);
imagesc(tfrsqticEXT*fs,ttc(1:basicTF.hop:end),log1p(abs(SSTxxLSE')/5e1)); colormap(1-gray);
ylim([tc(end)-0.6*twin/2 tc(end)]); %yticks([200 220 240]);
ylabel('Time (s)','rotation',270,'VerticalAlignment','top'); xlabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
set(gca,'yticklabelrotation',270) ;
set(gca,'xaxisLocation','top');

ax = [0.796 0.736] ;
ay = [0.724 0.724] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.276 0.216] ;  
ay = [0.670 0.670] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)

ax = [0.797 0.737] ;
ay = [0.103 0.103] ;
annotation('arrow',ax,ay,'Color','blue','linewidth',3)
ax = [0.281 0.221] ;  
ay = [0.069 0.069] ;
annotation('arrow',ax,ay,'Color','green','linewidth',3)

%STFT
figure;
subplot('Position',[0.055 0.12 0.435 0.865]);
imagesc(tc(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(Vx)/5e1)); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
subplot('Position',[0.56 0.12 0.435 0.865]);
imagesc(ttc(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(VxxLSE)/5e1)); xlim([0 tc(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;

dim = [.4625 0.12 .027 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;
dim = [.968 0.12 .027 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;

%% Performance index for SST and STFT

tSSTx = t(1:basicTF.hop:end);
tSSText = tt(1:basicTF.hop:end);
n0 = find(tSSText>=min(tSSTx),1,'first') ;
nf = find(tSSText<=max(tSSTx),1,'last') ;

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