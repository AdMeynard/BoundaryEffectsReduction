%% Display the results of InfluenceNoiseSine.m and ForecastAMFM_TBATS.r
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all;

load('../../Results/PerfNoise.mat')
L = size(VarXPm,2) ;
nbXP = length (Sigma) ;

%% Plot
ll = 1:L ;

for k = 1:nbXP
    sigma = Sigma(k) ;
    txtXP = ['Experimental \sigma = ',num2str(sigma,'%.0e')];
    txtTH = ['Theoretical    \sigma = ',num2str(sigma,'%.0e')];
    figure(1);
    plot(ll,MeanXPm(k,:)/sigma,'linewidth',2,'DisplayName',txtXP); grid on; hold on;
    figure(2);
    fig = semilogy(ll,VarXPm(k,:),'+','linewidth',2,'MarkerSize',8,'DisplayName',txtXP); grid on; hold on;
    semilogy(ll,VarTH(k,:),'-','color',fig.Color,'linewidth',2,'DisplayName',txtTH); hold on;
end
%legend('Experimental MSE[10]','Experimental MSE[100]','Experimemtal MSE[200]','Experimental MSE[500]','\sigma^2');
figure(1); xlabel('Forecasting sample index ($\ell$)','interpreter', 'latex','fontsize',28); ylabel('Normalized experimental forecasting bias ($\mu_{xp}[N-1+\ell]/sigma$)','interpreter','latex','fontsize',28);
legend show ; %ylim([-0.05 0.05]);
figure(2); xlabel('Forecasting sample index ($\ell$)','interpreter', 'latex','fontsize',28); ylabel('Forecasting variance ($\gamma[N-1+\ell,N-1+\ell]$)','interpreter','latex','fontsize',28);
lgd = legend('show') ;
lgd.FontSize = 18 ;
set(gca,'fontsize',24) ;
