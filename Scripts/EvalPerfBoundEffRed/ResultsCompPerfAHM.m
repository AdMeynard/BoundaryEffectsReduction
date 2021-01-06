%% Display the results of CompPerfAHM.m and ForecastAMFM_TBATS.r (Table I in the paper)
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

load('../../Results/PerfAHM') ;

dataTBATS = table2array( readtable('../../Results/PerfAHM_TBATS.csv','Range','B:D','TreatAsEmpty','NA') ) ;
BiasXP.TBATS = dataTBATS(:,1) ;
VarianceXP.TBATS = dataTBATS(:,2) ;
CPUtimeXP.TBATS = dataTBATS(~isnan(dataTBATS(:,3)),3) ;

MSE.LSE = BiasXP.LSE.^2 + VarianceXP.LSE ;
MSE.SYM = BiasXP.SYM.^2 + VarianceXP.SYM ;
MSE.EDMD = BiasXP.EDMD.^2 + VarianceXP.EDMD ;
MSE.GPR = BiasXP.GPR.^2 + VarianceXP.GPR ;
MSE.TBATS = BiasXP.TBATS.^2 + VarianceXP.TBATS ;

fprintf('======================= Forecasting Performance ===================\n')
fprintf(' ___________________________________________________________________\n')
fprintf('|    Extension   |  Computing  |                   MSE              |\n')
fprintf('|     Method     | time (sec.) |   Mean  |   SD  | t-test to SigExt |\n')
fprintf('|----------------|-------------|---------|-------|------------------|\n')
fprintf('|      SigExt    |     %.3f   |  %.3f  | %.3f |                  |\n', CPUtimeXP.LSE/2 , mean(MSE.LSE), std(MSE.LSE)) % time divided by 2 because 2 extensions are operated in CompPerfAHM
fprintf('| Symmetrization |     %.3f   |  %.3f  | %.3f |         %u        |\n', CPUtimeXP.SYM/2 , mean(MSE.SYM), std(MSE.SYM), ttest(MSE.SYM,MSE.LSE))
fprintf('|      EDMD      |    %.3f   |  %.3f  | %.3f |         %u        |\n', CPUtimeXP.EDMD/2 , mean(MSE.EDMD), std(MSE.EDMD), ttest(MSE.EDMD,MSE.LSE))
fprintf('|      GPR       |   %.3f   |  %.3f  | %.3f |         %u        |\n', CPUtimeXP.GPR/2 , mean(MSE.GPR), std(MSE.GPR), ttest(MSE.GPR,MSE.LSE))
fprintf('|      TBATS     |   %.3f   |  %.3f  | %.3f |         %u        |\n', CPUtimeXP.TBATS/2 , mean(MSE.TBATS), std(MSE.TBATS), ttest(MSE.TBATS,MSE.LSE))
fprintf('|________________|_____________|_________|_______|__________________|\n\n')