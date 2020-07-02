clear all; close all; clc;
load ../../Results/resultSucForPPGnoise ;

figure;
semilogy(SNR,NoiseForecastErr,'+','linewidth',2);
xlabel('SNR (dB)'); ylabel('Forecasting Variance');
grid on;
set(gca,'fontsize',24) ;

figure;
semilogy(SNR,NoiseTFR.STFT,'^',SNR,NoiseTFR.SST,'+',SNR,NoiseTFR.RS,'o','linewidth',2,'MarkerSize',8);
xlabel('SNR (dB)'); ylabel('Averaged OTD'); grid on;
legend({'STFT','SST','RS'},'location','southwest');
set(gca,'fontsize',24) ;