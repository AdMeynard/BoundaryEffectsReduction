clear all; close all; clc;
load ../../Results/resultSucForPPGnoise ;

figure;
loglog(StdNoise.^2,NoiseForecastErr,'+','linewidth',2);
xlabel('Noise Variance \sigma_{add}^2'); ylabel('Forecasting Variance');
grid on;
set(gca,'fontsize',24) ;

figure;
loglog(StdNoise.^2,NoiseTFR.SST,'+',StdNoise.^2,NoiseTFR.STFT,'^','linewidth',2,'MarkerSize',8); %,StdNoise.^2,NoiseRS,'+',
xlabel('Noise Variance \sigma_{add}^2'); ylabel('Averaged OTD'); grid on;
legend('SST','STFT');
set(gca,'fontsize',24) ;