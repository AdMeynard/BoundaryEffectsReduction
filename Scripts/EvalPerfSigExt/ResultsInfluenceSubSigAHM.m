%% Display the results of InfluenceSubsigLengthAHM.m
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;

load('../../Results/PerfSubsigLengthAHM') ;

BiasXP = abs(mean(MeanSigExt)) ;
MSEXP = VarSigExt ;
CPUtimeXP = mean(CPUtimeXP) ;

% figure;
% plot(extMval,BiasXP,'linewidth',2) ; axis tight ;
% xlabel('Subsignals length $M$','interpreter','latex') ;
% ylabel('Experimental Bias $\mu_{\mathrm{xp}}$','interpreter','latex') ;
% set(gca,'fontsize',24) ; grid on ;
% 
% figure;
% plot(extMval,MSEXP,'linewidth',2) ; axis tight ;
% xlabel('Subsignals length $M$','interpreter','latex') ;
% ylabel('Experimental MSE','interpreter','latex') ;
% set(gca,'fontsize',24) ; grid on ;


fprintf('============= Forecasting Performance =============\n')
fprintf(' ________________________________________________\n')
fprintf('|    Extension    |  Computing  |        MSE      |\n')
fprintf('|     Method      | time (sec.) |   Mean  |   SD  |\n')
fprintf('|-----------------|-------------|---------|-------|\n')

indM = 1 ;
for extM = extMval
    fprintf('| SigExt (M=%4i) |     %.3f   |  %.3f  | %.3f |\n', extM, CPUtimeXP(indM) , mean(MSEXP(:,indM)), std(MSEXP(:,indM)) ) ;
    indM = indM + 1 ;
end
indM = indM +1 ;

fprintf('|_________________|_____________|_________|_______|\n\n')