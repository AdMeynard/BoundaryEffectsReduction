clear all; close all; clc;
load ../../Results/resultSucForPPGnoise ;

figure;
loglog(StdNoise.^2,NoiseForecastErr,'+','linewidth',2);
xlabel('Noise Variance \sigma_{add}^2'); ylabel('Forecasting Variance');
grid on;
set(gca,'fontsize',24) ;

figure;
loglog(StdNoise.^2,NoiseTFR.SST,'+',StdNoise.^2,NoiseTFR.STFT,'^',StdNoise.^2,NoiseTFR.RS,'o','linewidth',2,'MarkerSize',8);
xlabel('Noise Variance \sigma_{add}^2'); ylabel('Averaged OTD'); grid on;
legend({'SST','RS','STFT'},'location','southeast');
set(gca,'fontsize',24) ;