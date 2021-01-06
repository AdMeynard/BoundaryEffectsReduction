%% Display the results of InfluenceSubsigLengthAHM.m
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

load('../../Results/PerfSubsigLengthAHM') ;

BiasXP = abs(mean(MeanSigExt)) ;
VarianceXP = mean(VarSigExt) ;

figure;
plot(extMval,BiasXP,'linewidth',2) ;
xlabel('Subsignals length $M$','interpreter','latex') ;
ylabel('Experimental Bias $\gamma_{\mathrm{xp}}$','interpreter','latex') ;
set(gca,'fontsize',24) ; grid on ;

figure;
plot(extMval,VarianceXP,'linewidth',2) ;
xlabel('Subsignals length $M$','interpreter','latex') ;
ylabel('Experimental Variance $\gamma_{\mathrm{xp}}$','interpreter','latex') ;
set(gca,'fontsize',24) ; grid on ;
