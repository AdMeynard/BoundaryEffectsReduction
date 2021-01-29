%% Run the extension methods, and the associated BoundEffRed TF representations on a EEG signal
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');

%% load PPG signal
addpath ../../Signals/
data = load('EEG');

xtot = data.EEG ;
fs = 200 ; % sampling frequency

n0 = 75e3;
nf = 85e3;
x = xtot(n0:nf);
mu = mean(x) ;
s = std(x) ;
x = (x - mu)/s ;

N = length(x);

%% Forecasting

HOP = 1 ;
extSEC = 1.5 ; % the extension is of extSEC second
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
fprintf('|       SigExt     |         %.2f          |\n', LSEtime)
fprintf('|  EDMD extension  |         %.2f          |\n', EDMDtime)
fprintf('|  GPR extension   |        %.2f          |\n', GPRtime)
fprintf('|__________________|_______________________|\n')

%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*round(extSEC*fs)) ;

xxTRUE = ( xtot((n0-L):(nf+L)) - mu ) / s ;
xxZP = [ zeros(L,1); x; zeros(L,1) ] ; % zero padding


figure;
plot(tt,xxLSE,tt,xxEDMD,tt,xxGPR,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
set(gca,'fontsize',14);
legend('LSE Estimated Extended signal','EDMD Estimated Extended signal','GPR Estimated Extended signal','Ground truth Extended signal','Original signal'); 
xlabel('Time (s)'); ylabel('Signals'); title('Time series'); axis tight;

%% SST
basicTF.hop = 10;
basicTF.win = 201;

% On the original signal
[~, SSTx, ~, tfrsqtic] = ConceFT_sqSTFT_C(double(x-mean(x)), 0, 0.15,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[~, SSTxxTRUE, ~, tfrsqticTRUE] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), 0, 0.15,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the zero-padded extended signal 
[~, SSTxxZP, ~, tfrsqticZP] = ConceFT_sqSTFT_C(double(xxZP-mean(xxZP)), 0, 0.15,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, SSTxxLSE, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), 0, 0.15,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, SSTxxEDMD, ~, ~] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), 0, 0.15,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, SSTxxGPR, ~, ~] = ConceFT_sqSTFT_C(double(xxGPR-mean(xxGPR)), 0, 0.15,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        
figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(SSTx)/1e1)); axis xy ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(SSTxxTRUE)/1e1)); xlim([0 t(end)]); axis xy ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxLSE)/1e1)); xlim([0 t(end)]); axis xy ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on LSE estimated long signal (short signal extended)');
subplot(2,2,4);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(SSTxxEDMD)/1e1)); xlim([0 t(end)]); axis xy ;
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on LSE estimated long signal (short signal extended)');

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

