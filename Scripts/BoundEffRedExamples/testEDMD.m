clear all; 
close all; clc

%% load THO signal
addpath ../../Signals/
data = load('THO.mat');
xtot = data.THO;
fs = 100 ; % sampling frequency

x = xtot(200e3:208e3) ;
m = mean(x) ;
s = std(x) ;
x = (x-m)/s ;
N = length(x);

%% Construct X and Y for Forecasting

HOP = 1 ;
extSEC = 12 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round(1.5*L) ; % dimension of embedding / signals length
extK = round( 2*extM ) ; % extP ;  % number of points to estimate A / size of datasets
if extK + extM >length(x) - 10 
    extK = extK/2 ; extM = extM/2 ;
end

X = [] ; Y = [] ;
for kk = 1: extK
    X(:,kk) = x(end-extK-extM+kk: HOP: end-extK+kk-1) ;
    Y(:,kk) = x(end-extK-extM+kk+1: HOP: end-extK+kk) ;
end

X = X.' ;
Y = Y.' ;

%% Koopman

sigma2 = 100 ;

tic;
[Xi,mu,phix] = approxKoopman(X,Y,sigma2) ;
toc;

%% Prediction / Extension

K = round(fs*extSEC);
Z = zeros(extK,K) ;

tmp = phix.';
for kk = 1:K
    tmp = mu.*tmp;
    Z(:,kk) =  tmp ; 
end
Z = real(Xi.' * Z);

xext = Z(end,:)' ;

xx = [x; xext];

xxTRUE = (xtot(200e3:(208e3+K)) -m)/s ;

t = linspace(0, (N-1)/fs, N) ;
tt = linspace(0, (N+K-1)/fs, N+K) ;
plot(tt,xxTRUE,tt,xx,t,x,'linewidth',2) ;
axis tight; grid on;
xlabel('Time (s)'); ylabel('Signal');
legend('Ground truth','Forecasting','Measured part')