%% Display the results of CompPerfAHM.m and ForecastAMFM_TBATS.r (Table I in the paper)
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

load('../../Results/PerfAHM_TBATS') ;
MSE.TBATS = VarianceXP.TBATS ;

load('../../Results/PerfAHM') ;

MSE.LSE =  VarianceXP.LSE ;
MSE.SYM = VarianceXP.SYM ;
MSE.EDMD = VarianceXP.EDMD ;
MSE.GPR = VarianceXP.GPR ;

fprintf('===================== Forecasting Performance ===================\n')
fprintf(' _______________________________________________________________\n')
fprintf('|    Extension   |  Computing  |              MSE               |\n')
fprintf('|     Method     | time (sec.) |  Mean |   SD  |t-test to SigExt|\n')
fprintf('|----------------|-------------|-------|-------|----------------|\n')

indM = 1 ;
for extM = extMval
    fprintf('|SigExt (M=%4i) |  %8.3f   | %5.3f | %5.3f |                |\n', extM, CPUtimeXP.LSE(indM) , mean(MSE.LSE(:,indM)), std(MSE.LSE(:,indM)) ) ;
    indM = indM + 1 ;
end
indM = indM +1 ;

MSEref = MSE.LSE(:,2) ;

fprintf('| Symmetrization |  %8.3f   | %5.3f | %5.3f |        %u       |\n', CPUtimeXP.SYM , mean(MSE.SYM), std(MSE.SYM), ttest(MSE.SYM,MSEref))
fprintf('|      EDMD      |  %8.3f   | %5.3f | %5.3f |        %u       |\n', CPUtimeXP.EDMD , mean(MSE.EDMD), std(MSE.EDMD), ttest(MSE.EDMD,MSEref))
fprintf('|      GPR       |  %8.3f   | %5.3f | %5.3f |        %u       |\n', CPUtimeXP.GPR , mean(MSE.GPR), std(MSE.GPR), ttest(MSE.GPR,MSEref))

load('../../Results/PerfAHM_TBATS') ;
MSE.TBATS = VarianceXP.TBATS ;

fprintf('|      TBATS     |  %8.3f   | %5.3f | %5.3f |                |\n', CPUtimeXP.TBATS , mean(MSE.TBATS), std(MSE.TBATS))
fprintf('|________________|_____________|_______|_______|________________|\n\n')

load('../../Results/PerfAHM') ;
OTDref = OTD.STFT.S ;

fprintf('========= BoundEffRed =============\n')
fprintf(' __________________________________\n')
fprintf('|    Extension   |Performance Index|\n')
fprintf('|     Method     |  Mean  |   SD   |\n')
fprintf('|----------------|--------|--------|\n')

indM = 1 ;
for extM = extMval
    fprintf('|SigExt (M=%4i) | %6.4f | %6.4f |\n', extM, mean(OTD.STFT.LSE(:,indM)./OTDref'), std(OTD.STFT.LSE(:,indM)./OTDref')) ;
    indM = indM + 1 ;
end
indM = indM +1 ;

fprintf('| Symmetrization | %6.4f | %6.4f |\n', mean(OTD.STFT.SYM./OTDref), std(OTD.STFT.SYM./OTDref))
fprintf('|      EDMD      | %6.4f | %6.4f |\n', mean(OTD.STFT.EDMD./OTDref), std(OTD.STFT.EDMD./OTDref))
fprintf('|      GPR       | %6.4f | %6.4f |\n', mean(OTD.STFT.GPR./OTDref), std(OTD.STFT.GPR./OTDref))

load('../../Results/PerfAHM_TBATS') ;
fprintf('|      TBATS     | %6.4f | %6.4f |\n', mean(OTD.STFT.TBATS./OTDref), std(OTD.STFT.TBATS./OTDref))
fprintf('|________________|________|________|\n\n')