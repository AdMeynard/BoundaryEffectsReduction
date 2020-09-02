%% Display the results of InfluenceNoiseSine.m and InfluenceSizeDatasetSine.m
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all;

%% Influence of the noise level
load('../../Results/PerfNoise.mat')

figure;
% subplot(2,1,1) ;
loglog(Sigma.^2,abs(BiasXPm(:,1)),Sigma.^2,abs(BiasXPm(:,10)),Sigma.^2,abs(BiasXPm(:,100)),'linewidth',2);
axis tight; ylim([1e-10 1e3]); grid on;
xlabel('Noise Variance $\sigma^2$','interpreter','latex') ; ylabel('Experimental Bias $\mu_{\mathrm{xp}}$','interpreter','latex') ;
legend({'$\ell=1$','$\ell=10$','$\ell=100$'},'interpreter','latex') ;
set(gca,'fontsize',24) ;

figure;
% subplot(2,1,2) ;
loglog(Sigma.^2,VarXPm(:,1),Sigma.^2,VarXPm(:,10),Sigma.^2,VarXPm(:,100),'linewidth',2);
axis tight; ylim([1e-15 1e3]); grid on;
xlabel('Noise Variance $\sigma^2$','interpreter','latex') ; ylabel('Experimental Variance $\gamma_{\mathrm{xp}}$','interpreter','latex') ;
legend({'$\ell=1$','$\ell=10$','$\ell=100$'},'interpreter','latex') ;
set(gca,'fontsize',24) ;

%% Influence of the dataset size
load('../../Results/PerfSizeDataset.mat')

figure;
% subplot(2,1,1) ;
loglog(KK,abs(BiasXPm(:,1)),KK,abs(BiasXPm(:,10)),KK,abs(BiasXPm(:,100)),'linewidth',2);
axis tight; ylim([1e-8 1e0]); grid on;
xlabel('Dataset Size $K$','interpreter','latex') ; ylabel('Experimental Bias $\mu_{\mathrm{xp}}$','interpreter','latex') ;
legend({'$\ell=1$','$\ell=10$','$\ell=100$'},'interpreter','latex') ; xticks([300 500 1000 2000]);
set(gca,'fontsize',24) ;

figure;
% subplot(2,1,2) ;
loglog(KK,VarXPm(:,1),KK,VarXPm(:,10),KK,VarXPm(:,100),'linewidth',2);
axis tight; grid on;
xlabel('Dataset Size $K$','interpreter','latex') ; ylabel('Experimental Variance $\gamma_{\mathrm{xp}}$','interpreter','latex') ;
legend({'$\ell=1$','$\ell=10$','$\ell=100$'},'interpreter','latex') ; xticks([300 500 1000 2000]);
set(gca,'fontsize',24) ;
