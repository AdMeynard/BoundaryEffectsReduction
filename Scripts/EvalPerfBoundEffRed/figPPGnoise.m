%% Display the results of EvalPPGnoise.m
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;
load ../../Results/resultSucForPPGnoise ;

%% SNR
figure;
semilogy(SNR,NoiseForecastErr,'+','linewidth',2);
xlabel('SNR (dB)'); ylabel('Forecasting Variance');
grid on;
set(gca,'fontsize',24) ;

figure;
plot(SNR,NoiseTFR.STFT,'^',SNR,NoiseTFR.SST,'+',SNR,NoiseTFR.RS,'o','linewidth',2,'MarkerSize',10);
ylim([0 0.6]) ;
xlabel('SNR (dB)','interpreter','latex'); ylabel('Averaged performance index $D$','interpreter','latex'); grid on;
legend('STFT','SST','RS');
set(gca,'fontsize',24) ;