clear all; close all;

load('../../Results/MSEnoise.mat')
extM = 750 ;
nbXP = length(MSE100m);
Sigma = linspace(5e-3,1e-1,nbXP) ;

%% Plot
tmp = 2/extM ;
MSE200th = Sigma.^2*(1+tmp*(1+tmp)^(2*200-2)) ;

figure;
plot(Sigma.^2,MSE10m,Sigma.^2,MSE100m,Sigma.^2,MSE200m,Sigma.^2,MSE500m,Sigma.^2, Sigma.^2,'k--',Sigma.^2, MSE200th,'k--','linewidth',2); grid on;
legend('Experimemtal MSE[10]','Experimemtal MSE[100]','Experimemtal MSE[200]','Experimemtal MSE[500]','\sigma^2'); 
xlabel('\sigma^2'); ylabel('Experimemtal MSE'); title('MSE in function of the noise level');