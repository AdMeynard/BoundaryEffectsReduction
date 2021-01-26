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

%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(0, (N-1+L)/fs, N+L) ;

xxTRUE = ( xtot(n0:(nf+L)) - mu ) / s ;

figure;
plot(tt,xxLSE,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
ylabel('Signals'); axis tight; xlim([t(end)-1.5*extSEC tt(end)]);
xticks([]);
legend({'{\sf SigExt} extension','Ground truth extension','Original signal'},'location','northwest','interpreter','latex');
set(gca,'fontsize',18);

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
% On the estimated extended signal 
[VxxLSE, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE, fmin, fmax,df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

figure;
subplot('Position',[0.055 0.12 0.435 0.865]);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(SSTx)/5e1)); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
subplot('Position',[0.56 0.12 0.435 0.865]);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;

dim = [.4325 0.12 .0575 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;
dim = [.938 0.12 .0575 0.865] ;
annotation('rectangle',dim,'Color','red','linewidth',3,'linestyle','--') ;

figure;
subplot('Position',[0.062 0.516 0.9 0.37]);
imagesc(tfrsqtic*fs,t(1:basicTF.hop:end),log1p(abs(SSTx')/5e1)); colormap(1-gray);
ylim([t(end)-0.6*twin/2 t(end)]); %yticks([200 220 240]);
ylabel('Time (s)','rotation',270,'VerticalAlignment','top'); xlabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
set(gca,'yticklabelrotation',270) ;
set(gca,'xaxisLocation','top')
subplot('Position',[0.062 0.02 0.9 0.37]);
imagesc(tfrsqticEXT*fs,tt(1:basicTF.hop:end),log1p(abs(SSTxxLSE')/5e1)); colormap(1-gray);
ylim([t(end)-0.6*twin/2 t(end)]); %yticks([200 220 240]);
ylabel('Time (s)','rotation',270,'VerticalAlignment','top'); xlabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
set(gca,'yticklabelrotation',270) ;
set(gca,'xaxisLocation','top');

figure;
subplot('Position',[0.055 0.12 0.435 0.865]);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(Vx)/5e1)); axis xy ; colormap(1-gray);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',24) ;
subplot('Position',[0.56 0.12 0.435 0.865]);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(VxxLSE)/5e1)); xlim([0 t(end)]); axis xy ; colormap(1-gray);
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