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
%legend('Experimemtal MSE[10]','Experimemtal MSE[100]','Experimemtal MSE[200]','Experimemtal MSE[500]','\sigma^2');
figure(1); xlabel('Forecasting interval'); ylabel('Experimental Forecasting Mean');

BndVarTH = 2*tmp*(1+tmp).^(2*ll-2) ;
figure(2); 
semilogy(ll,BndVarTH,'--','linewidth',2);
xlabel('Forecasting interval'); ylabel('Relative Experimental Forecasting Variance (Var/\sigma^2)');