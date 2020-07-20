%% Display the results of InfluenceSizeDatasetSine.m
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all;

load('../../Results/PerfSizeDataset.mat')
L = size(VarXPm,2) ;
nbXP = length (KK) ;

%% Plot
ll = 1:L ;

for k = 1:nbXP % k must be sufficiently large to observe the asymptotic behavior
    txtXP = ['Experimental K = ',num2str(KK(k),'%.0f')];
    txtTH = ['Theoretical    K = ',num2str(KK(k),'%.0f')];
    figure(1);
    plot(ll,MeanXPm(k,:),'linewidth',2,'DisplayName',txtXP); grid on; hold on;
    figure(2);
    fig = semilogy(ll,VarXPm(k,:),'+','linewidth',2,'DisplayName',txtXP); grid on; hold on;
    semilogy(ll,VarTH(k,:),'-','color',fig.Color,'linewidth',2,'DisplayName',txtTH); hold on;
end
figure(1); xlabel('Forecasting sample index ($\ell$)','interpreter', 'latex','fontsize',28); ylabel('Experimental forecasting bias ($\mu_{xp}[\ell]$)','interpreter','latex','fontsize',28);
ylim([-0.01 0.01]); legend show ;
figure(2); xlabel('Forecasting sample index ($\ell$)','interpreter', 'latex','fontsize',28); ylabel('Forecasting variance','interpreter','latex','fontsize',28);
lgd = legend('show') ;
lgd.FontSize = 18 ;
set(gca,'fontsize',24) ;
