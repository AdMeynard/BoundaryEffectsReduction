clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');
addpath('../../Signals/');

%% load THO signal
data = load('THO.mat');
xtot = data.THO;
fs = 100 ; % sampling frequency

N = 6e3; % subsignals length

% k = 9 ;
k = 213 ;
n0 = k*N+1 ;
n1 = (k+1)*N ;
x = xtot(n0:n1) ;
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

N = length(x);

%% Forecasting

HOP = 1 ;
extSEC = 7 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round(1.5*L) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets
if extK + extM >length(x) - 10 
    extK = extK/2 ; extM = extM/2 ;
end

method.name = 'lseV' ;
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

method.name = 'symmetric' ;
tic; 
xxSYM = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
SYMtime = toc;

%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*round(extSEC*fs)) ;

xxTRUE = ( xtot((n0-L):(n1+L)) - mu ) / s ;
xxZP = [ zeros(L,1); x; zeros(L,1) ] ; % zero padding

xxLSE = flipud(xxLSE);
xxEDMD = flipud(xxEDMD);
xxGPR = flipud(xxGPR);
xxSYM = flipud(xxSYM);
xxTRUE = flipud(xxTRUE);
x = flipud(x) ;

figure;
subplot(2,2,[1 3]);
plot(tt,xxLSE,'Color','#0072BD','linewidth',2); hold on;
plot(tt,xxEDMD,'Color','#7E2F8E','linewidth',2)
plot(tt,xxSYM,'Color','#77AC30','linewidth',2)
plot(tt,xxTRUE,'Color','#D95319','linewidth',2,'linestyle','--')
plot(t,x,'Color','#EDB120','linewidth',2)
grid on;
legend({'{\sf SigExt} extension','Symmetric extension','EDMD extension','Ground truth extension','Original signal'},'interpreter','latex','location','northwest'); %,'GPR extension' 
xlabel('Time (s)'); ylabel('Signals'); axis tight; xlim([tt(end)-25 tt(end)]); %ylim([-1.5 1.5]);
set(gca,'FontSize',20); %xticklabels([]) ;


%% SST
basicTF.fmin = 0/fs ;
basicTF.fmax = 1/fs ;
basicTF.df = 1e-5 ;
basicTF.hop = 10;
basicTF.win = 2*L+1;


[~, SSTx, ~, tfrsqticX] = ConceFT_sqSTFT_C(x, basicTF.fmin, basicTF.fmax,...
            basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        
[~, SSTxxLSE, ~, tfrsqticXX] = ConceFT_sqSTFT_C(xxLSE, basicTF.fmin, basicTF.fmax,...
    basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

tSST = t(1:basicTF.hop:end) ;
ttSST = tt(1:basicTF.hop:end) ;

subplot(2,2,2);
imagesc(tSST,tfrsqticX*fs,log1p(abs(SSTx)/5e1)); colormap(1-gray); axis xy ;
xlim([tSST(end)-15 tSST(end)])
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'FontSize',20);
subplot(2,2,4);
imagesc(ttSST,tfrsqticXX*fs,log1p(abs(SSTxxLSE)/5e1)); colormap(1-gray); axis xy ;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
xlim([tSST(end)-15 tSST(end)])
set(gca,'FontSize',20);