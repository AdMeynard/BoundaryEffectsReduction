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

fprintf('======== Forecasting Performance ========\n')
fprintf(' ________________________________________\n')
fprintf('|    Extension   |  Computing  | Averaged|\n')
fprintf('|     Method     | time (sec.) |   MSE   |\n')
fprintf('|----------------|-------------|---------|\n')
fprintf('|      LSE       |     %.3f   |  %.3f  |\n', CPUtimeXP.LSE/2 , mean(MSE.LSE)) % time diveded by 2 because 2 extensions are operated in CompPerfAHM
fprintf('| Symmetrization |     %.3f   |  %.3f  |\n', CPUtimeXP.SYM/2 , mean(MSE.SYM))
fprintf('|      EDMD      |    %.3f   |  %.3f  |\n', CPUtimeXP.EDMD/2 , mean(MSE.EDMD))
fprintf('|      GPR       |   %.3f   |  %.3f  |\n', CPUtimeXP.GPR/2 , mean(MSE.GPR))
fprintf('|      TBATS     |   %.3f   |  %.3f  |\n', CPUtimeXP.TBATS/2 , mean(MSE.TBATS))
fprintf('|________________|_____________|_________|\n\n')