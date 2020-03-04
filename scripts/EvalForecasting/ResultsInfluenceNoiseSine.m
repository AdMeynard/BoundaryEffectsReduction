clear all; close all;

load('../../Results/MSEnoise.mat')

nbXP = length(MSE100m);
Sigma = linspace(5e-3,1e-1,nbXP) ;

%% Plot

figure;
plot(Sigma.^2,MSE10m,Sigma.^2,MSE100m,Sigma.^2,MSE200m,Sigma.^2,MSE500m,Sigma.^2,2 *Sigma.^2,'--','linewidth',2); grid on;
legend('Experimemtal MSE[10]','Experimemtal MSE[100]','Experimemtal MSE[200]','Experimemtal MSE[500]','2\sigma^2'); 
xlabel('\sigma^2'); ylabel('Experimemtal MSE'); title('MSE in function of the noise level');