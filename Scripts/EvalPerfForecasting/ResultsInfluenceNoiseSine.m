clear all; close all;

load('../../Results/PerfNoise.mat')
L = size(VarXPm,2) ;
nbXP = length (Sigma) ;

%% Plot
tmp = 2/extM ;
ll = 1:L ;


for k = 1:nbXP
    sigma = Sigma(k) ;
    figure(1);
    plot(ll,MeanXPm(k,:),'linewidth',2); grid on; hold on;
    figure(2);
    semilogy(ll,VarXPm(k,:)/sigma^2,'linewidth',2); grid on; hold on; %,ll,BndVarTH,'--',
end
%legend('Experimental MSE[10]','Experimental MSE[100]','Experimemtal MSE[200]','Experimental MSE[500]','\sigma^2');
figure(1); xlabel('Forecasting sample index ($\ell$)','interpreter', 'latex','fontsize',15); ylabel('Experimental forecasting bias ($\mu_{xp}[\ell]$)','interpreter','latex','fontsize',15);
ylim([-0.05 0.05]);
figure(2); xlabel('Forecasting sample index ($\ell$)','interpreter', 'latex','fontsize',15); ylabel('Normalized experimental forecasting variance ($\gamma_{xp}[\ell]/\sigma^2$)','interpreter','latex','fontsize',15);
ylim([0.2 1]);
% BndVarTH = 2*tmp*(1+tmp).^(2*ll-2) ;
% semilogy(ll,BndVarTH,'--','linewidth',2);
