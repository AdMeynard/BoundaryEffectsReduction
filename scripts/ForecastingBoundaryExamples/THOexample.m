clear all; close all; clc;
addpath(genpath('../../SST'));
addpath('../../algorithm/');
addpath('../../Signals/');

%% load THO signal
data = load('THO.mat');
xtot = data.THO;
fs = 100 ; % sampling frequency

x = xtot(204e3:208e3);
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

%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*round(extSEC*fs)) ;

xxTRUE = ( xtot((204e3-L):(208e3+L)) - mu ) / s ;
xxZP = [ zeros(L,1); x; zeros(L,1) ] ; % zero padding

figure;
plot(tt,xxLSE,tt,xxEDMD,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
legend('LSE Estimated Extended signal','EDMD Estimated Extended signal','Ground truth Extended signal','Original signal'); 
xlabel('Time (s)'); ylabel('Signals'); title('Time series');

%% SST
basicTF.hop = 10;
basicTF.win = 1501;

% On the original signal
[~, ~, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(double(x-mean(x)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[~, ~, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the zero-padded extended signal 
[~, ~, SSTxxZP, ~, tfrsqticZP] = ConceFT_sqSTFT_C(double(xxZP-mean(xxZP)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, ~, SSTxxLSE, ~, ~] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, ~, SSTxxEDMD, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), 0, 0.01,...
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
%% Optimal transport distance
tmp = tt(1:basicTF.hop:end) ;
tmp = (tmp>=0) & (tmp<=t(end)) ;

tfrsqTRUEw = SSTxxTRUE(:,tmp) ;
tfrsqZPw = SSTxxZP(:,tmp) ;
tfrsqLSEw = SSTxxLSE(:,tmp) ;
tfrsqEDMDw = SSTxxEDMD(:,tmp) ;

otdZP = slicedOT(tfrsqZPw, tfrsqTRUEw) ;
otdLSE = slicedOT(tfrsqLSEw, tfrsqTRUEw) ;
otdEDMD = slicedOT(tfrsqEDMDw, tfrsqTRUEw) ;

fprintf(' ______________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) |    OTD    |\n')
fprintf('|======================================================|\n')
fprintf('|  Zero-padding    |         %.2f          | %.3e |\n', 0, otdZP)
fprintf('|  LSE extension   |         %.2f          | %.3e |\n', LSEtime, otdLSE)
fprintf('|  EDMD extension  |         %.2f          | %.3e |\n', EDMDtime, otdEDMD)
fprintf('|__________________|_______________________|___________|\n')

