clear all; close all; clc;
addpath(genpath('../SST'));
addpath('../algorithm/');

%% load PPG signal
addpath ../Signals/
data = load('PPGsig');

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


xx = SigExtension(x,fs,HOP,extK,extM,extSEC,'lseV');

%% Plot
t = linspace(0, (N-1)/fs, N) ;
tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*round(extSEC*fs)) ;

xxTRUE = ( xtot((30e3-L):(34e3+L)) - mu ) / s ;
xxZP = [ zeros(L,1); x; zeros(L,1) ] ; % zero padding


figure;
plot(tt,xx,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
set(gca,'fontsize',14);
legend('Estimated Extended signal','Ground truth Extended signal','Original signal'); 
xlabel('Time (s)'); ylabel('Signals'); title('Time series'); axis tight;

%% SST
basicTF.hop = 10;
basicTF.win = 1501;

% On the original signal
[~, ~, tfrsq3, ConceFT3, tfrsqtic] = ConceFT_sqSTFT_C(double(x-mean(x)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the true extended signal 
[~, ~, tfrsq3TRUE, ConceFT3TRUE, tfrsqticTRUE] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the zero-padded extended signal 
[~, ~, tfrsq3ZP, ConceFT3ZP, tfrsqticZP] = ConceFT_sqSTFT_C(double(xxZP-mean(xxZP)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
% On the estimated extended signal 
[~, ~, tfrsq3EXT, ConceFT3EXT, tfrsqticEXT] = ConceFT_sqSTFT_C(double(xx-mean(xx)), 0, 0.01,...
            1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

figure;
subplot(2,2,1);
imagesc(t(1:basicTF.hop:end),tfrsqtic*fs,log1p(abs(tfrsq3)/1e2));
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original short signal');
subplot(2,2,2);
imagesc(tt(1:basicTF.hop:end),tfrsqticTRUE*fs,log1p(abs(tfrsq3TRUE)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on original long signal');
subplot(2,2,3);
imagesc(tt(1:basicTF.hop:end),tfrsqticZP*fs,log1p(abs(tfrsq3ZP)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on zero-padded long signal (short signal extended)');
subplot(2,2,4);
imagesc(tt(1:basicTF.hop:end),tfrsqticEXT*fs,log1p(abs(tfrsq3EXT)/1e2)); xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('SST on estimated long signal (short signal extended)');

%save('results','tfrsq3','tfrsq3EXT');
%% Optimal transport distance
tmp = tt(1:basicTF.hop:end) ;
tmp = (tmp>=0) & (tmp<=t(end)) ;

tfrsq3TRUEw = tfrsq3TRUE(:,tmp) ;
tfrsq3ZPw = tfrsq3ZP(:,tmp) ;
tfrsq3EXTw = tfrsq3EXT(:,tmp) ;

OTDshort = slicedOT(tfrsq3(:,1:(end-1)), tfrsq3TRUEw) ;
OTDZP = slicedOT(tfrsq3ZPw, tfrsq3TRUEw) ;
OTDEXT = slicedOT(tfrsq3EXTw, tfrsq3TRUEw) ;

fprintf(' Extension Method |    OTD    | \n')
fprintf('  Short signal    | %.3e |\n', OTDshort)
fprintf('  Zero-padding    | %.3e |\n', OTDZP)
fprintf('  LSE extension   | %.3e |\n', OTDEXT)

