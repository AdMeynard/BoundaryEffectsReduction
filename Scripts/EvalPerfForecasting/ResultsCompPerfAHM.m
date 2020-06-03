load('../../Results/PerfAHM') ;

dataTBATS = table2array( readtable('../../Results/PerfAHM_TBATS.csv','Range','B:D','TreatAsEmpty','NA') ) ;
BiasXP.TBATS = dataTBATS(:,1) ;
VarianceXP.TBATS = dataTBATS(:,2) ;
CPUtimeXP.TBATS = dataTBATS(~isnan(dataTBATS(:,3)),3) ;

MSE.LSE = BiasXP.LSE.^2 + VarianceXP.LSE ;
MSE.EDMD = BiasXP.EDMD.^2 + VarianceXP.EDMD ;
MSE.GPR = BiasXP.GPR.^2 + VarianceXP.GPR ;

fprintf('===================Forecasting Performance=========\n')
fprintf(' _________________________________________________\n')
fprintf('| Extension |  Computing  |   Forecasting Error   |\n')
fprintf('|  Method   | time (sec.) |    Mean   |  Variance |\n')
fprintf('|-----------|-------------|-----------|-----------|\n')
fprintf('|   LSE     |     %.2f    | %.3e | %.3e |\n', CPUtimeXP.LSE , median(MSE.LSE), var(MSE.LSE))
fprintf('|   EDMD    |     %.2f    | %.3e | %.3e |\n', CPUtimeXP.EDMD , median(MSE.EDMD), var(MSE.EDMD))
fprintf('|   GPR     |   %.2f    | %.3e | %.3e |\n', CPUtimeXP.GPR , median(MSE.GPR), var(MSE.GPR))
fprintf('|   TBATS   |   %.2f    | %.3e | %.3e |\n', CPUtimeXP.TBATS , median(MSE.TBATS), var(MSE.TBATS))
fprintf('|___________|_____________|___________|___________|\n\n')