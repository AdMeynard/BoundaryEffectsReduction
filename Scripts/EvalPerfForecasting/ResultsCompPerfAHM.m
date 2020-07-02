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

fprintf('===================Forecasting Performance=========\n')
fprintf(' ______________________________________________________\n')
fprintf('|    Extension   |  Computing  |   Forecasting Error   |\n')
fprintf('|     Method     | time (sec.) |    Mean   |  Variance |\n')
fprintf('|----------------|-------------|-----------|-----------|\n')
fprintf('|      LSE       |     %.3f   | %.3e | %.3e |\n', CPUtimeXP.LSE , mean(MSE.LSE), var(MSE.LSE))
fprintf('| Symmetrization |     %.3f   | %.3e | %.3e |\n', CPUtimeXP.SYM , mean(MSE.SYM), var(MSE.SYM))
fprintf('|      EDMD      |    %.3f   | %.3e | %.3e |\n', CPUtimeXP.EDMD , mean(MSE.EDMD), var(MSE.EDMD))
fprintf('|      GPR       |   %.3f   | %.3e | %.3e |\n', CPUtimeXP.GPR , mean(MSE.GPR), var(MSE.GPR))
fprintf('|      TBATS     |  %.3f   | %.3e | %.3e |\n', CPUtimeXP.TBATS , mean(MSE.TBATS), var(MSE.TBATS))
fprintf('|________________|_____________|___________|___________|\n\n')