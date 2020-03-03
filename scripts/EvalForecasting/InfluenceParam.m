clear all;
close all; clc;
addpath('../algorithm/');
addpath(genpath('../SST'));

%% load signal
SigTyp = 'PPG' ;

switch SigTyp
    case 'THO'
        addpath ../NonStationaritySleepDataSet/1
        data = load('THO.mat');
        xtot = data.THO;
        fs = 100 ; % sampling frequency
        
        x = xtot(204e3:208e3);
        mu = mean(x) ;
        s = std(x) ;
        x = (x - mu)/s ;

        extSEC = 2 ; % the extension is of extSEC second
        L = round(extSEC*fs) ;
        xxTRUE = ( xtot((204e3-L):(208e3+L)) - mu ) / s ;

    case 'PPG'
        addpath ../Signals/
        data = load('PPGsig');

        xtot = data.PPG;
        fs = data.fs; % sampling frequency

        x = xtot(30e3:34e3);
        mu = mean(x) ;
        s = std(x) ;
        x = (x - mu)/s ;
        
        extSEC = 2 ; % the extension is of extSEC second
        L = round(extSEC*fs) ;
        xxTRUE = ( xtot((30e3-L):(34e3+L)) - mu ) / s ;
end

N = length(x);
t = linspace(0, (N-1)/fs, N);

figure;
plot(t,x,'linewidth',2);
set(gca,'fontsize',16);
xlabel('Time (s)'); ylabel('PPG signal'); grid on; axis tight;

%% Influence of M (Fix K) 
HOP = 1 ;
Mh = round( linspace(2,5*L,100) );
extK = 10*L ;

m = 1;
for extM = Mh
    xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,'lse');
    xxLSEV = SigExtension(x,fs,HOP,extK,extM,extSEC,'lseV');
%     xxDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,'dmd');
    forecastErrLSE_m(m) = sqrt( (1/(2*L)) * sum(abs(xxLSE - xxTRUE).^2) ) ;
    forecastErrLSEV_m(m) = sqrt( (1/(2*L)) * sum(abs(xxLSEV - xxTRUE).^2) ) ;
%     forecastErrDMD_m(m) = sqrt((1/(2*L)) * sum(abs(xxDMD - xxTRUE).^2) ) ;
    m = m+1 ;
end

figure;
plot( Mh/L,forecastErrLSE_m, 'linewidth',2); % , Mh/L,forecastErrLSEV_m
set(gca,'fontsize',16);
%legend( 'LSE', 'LSE Vector');
xlabel('$M/L$', 'interpreter','latex'); ylabel('MSE'); grid on;

%% Influence of K (Fix M) 

extM = round( 1.5*L );
Kh = round( linspace(1.05,5,100)*extM );

k = 1;
for extK = Kh
    xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,'lse');
    xxLSEV = SigExtension(x,fs,HOP,extK,extM,extSEC,'lseV');
    %xxDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,'dmd');
    forecastErrLSE_k(k) = sqrt((1/(2*L)) * sum(abs(xxLSE - xxTRUE).^2) ) ;
    forecastErrLSEV_k(k) = sqrt( (1/(2*L)) * sum(abs(xxLSEV - xxTRUE).^2) ) ;
    %forecastErrDMD_k(k) = sqrt((1/(2*L)) * sum(abs(xxDMD - xxTRUE).^2) ) ;
    k = k+1 ;
end


figure;
plot( Kh/L,forecastErrLSE_k, 'linewidth',2 ); %, Kh/L,forecastErrLSEV_k
set(gca,'fontsize',16);
%legend( 'LSE', 'LSE Vector');
xlabel('$K/L$', 'interpreter','latex'); ylabel('MSE'); grid on;
xlim([0 Kh(end)/L]); ylim([0 max(forecastErrLSE_k)]);
%% Influence of additive GWN

extM = round( 1.5*L );
extK = round( 2.5*extM );

sigmah = linspace(0,1,100) ;

n = 1;
for sigma = sigmah
    noise = sigma*randn(N,1) ;
    xn = x + noise;
    
    xxLSE = SigExtension(xn,fs,HOP,extK,extM,extSEC,'lse');
    xxLSEV = SigExtension(xn,fs,HOP,extK,extM,extSEC,'lseV');
%     xxDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,'dmd');
    forecastErrLSE_n(n) = sqrt( (1/(2*L)) * sum(abs(xxLSE - xxTRUE).^2) ) ;
    forecastErrLSEV_n(n) = sqrt((1/(2*L)) * sum(abs(xxLSEV - xxTRUE).^2) ) ;
%     forecastErrDMD_n) = sqrt((1/(2*L)) * sum(abs(xxDMD - xxTRUE).^2) ) ;
    n = n+1 ;
end

figure;
plot( sigmah.^2,forecastErrLSE_n.^2 - sigmah.^2 , 'linewidth',2 ); %, sigmah,sqrt( forecastErrLSEV_n.^2 - sigmah.^2 ),
set(gca,'fontsize',16); grid on;
%legend( 'LSE', 'LSE Vector');
xlabel('$\sigma^2$','interpreter','latex'); ylabel('$MSE - \sigma^2$','interpreter','latex');
