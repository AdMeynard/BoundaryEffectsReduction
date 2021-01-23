%% Display the results of CompPerfAHM.m and ForecastAMFM_TBATS.r (Table I in the paper)
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

load('../../Results/PerfAHM') ;

% dataTBATS = table2array( readtable('../../Results/PerfAHM_TBATS.csv','Range','B:D','TreatAsEmpty','NA') ) ;
% BiasXP.TBATS = dataTBATS(:,1) ;
% MSE.TBATS = dataTBATS(:,2) ;
% CPUtimeXP.TBATS = dataTBATS(~isnan(dataTBATS(:,3)),3) ;

MSE.LSE =  VarianceXP.LSE ;
MSE.SYM = VarianceXP.SYM ;
MSE.EDMD = VarianceXP.EDMD ;
MSE.GPR = VarianceXP.GPR ;
% MSE.TBATS = dataTBATS(:,1).^2 + dataTBATS(:,2); %VarianceXP.TBATS ;

fprintf('===================== Forecasting Performance ===================\n')
fprintf(' _______________________________________________________________\n')
fprintf('|    Extension   |  Computing  |              MSE               |\n')
fprintf('|     Method     | time (sec.) |  Mean |   SD  |t-test to SigExt|\n')
fprintf('|----------------|-------------|-------|-------|----------------|\n')

indM = 1 ;
for extM = extMval
    fprintf('|SigExt (M=%4i) |   %7.3f   | %5.3f | %5.3f |                |\n', extM, CPUtimeXP.LSE(indM) , mean(MSE.LSE(:,indM)), std(MSE.LSE(:,indM)) ) ;
    indM = indM + 1 ;
end
indM = indM +1 ;

fprintf('| Symmetrization |   %7.3f   | %5.3f | %5.3f |        %u       |\n', CPUtimeXP.SYM , mean(MSE.SYM), std(MSE.SYM), ttest(MSE.SYM,MSE.LSE(:,2)))
fprintf('|      EDMD      |   %7.3f   | %5.3f | %5.3f |        %u       |\n', CPUtimeXP.EDMD , mean(MSE.EDMD), std(MSE.EDMD), ttest(MSE.EDMD,MSE.LSE(:,2)))
fprintf('|      GPR       |   %7.3f   | %5.3f | %5.3f |        %u       |\n', CPUtimeXP.GPR , mean(MSE.GPR), std(MSE.GPR), ttest(MSE.GPR,MSE.LSE(:,2)))
% fprintf('|      TBATS     |   %7.3f   | %5.3f | %5.3f |         %u        |\n', CPUtimeXP.TBATS , mean(MSE.TBATS), std(MSE.TBATS), ttest(MSE.TBATS,MSE.LSE))
fprintf('|________________|_____________|_______|_______|________________|\n\n')

fprintf('========= BoundEffRed =============\n')
fprintf(' __________________________________\n')
fprintf('|    Extension   |Performance Index|\n')
fprintf('|     Method     |  Mean  |   SD   |\n')
fprintf('|----------------|--------|--------|\n')

indM = 1 ;
for extM = extMval
    fprintf('|SigExt (M=%4i) | %6.4f | %6.4f |\n', extM, mean(OTD.STFT.LSE(:,indM)./OTD.STFT.S'), std(OTD.STFT.LSE(:,indM)./OTD.STFT.S')) ;
    indM = indM + 1 ;
end
indM = indM +1 ;

fprintf('| Symmetrization | %6.4f | %6.4f |\n', mean(OTD.STFT.SYM./OTD.STFT.S), std(OTD.STFT.SYM./OTD.STFT.S))
fprintf('|      EDMD      | %6.4f | %6.4f |\n', mean(OTD.STFT.EDMD./OTD.STFT.S), std(OTD.STFT.EDMD./OTD.STFT.S))
fprintf('|      GPR       | %6.4f | %6.4f |\n', mean(OTD.STFT.GPR./OTD.STFT.S), std(OTD.STFT.GPR./OTD.STFT.S))
% fprintf('|      TBATS     | %6.4f | %6.4f |\n', mean(OTD.STFT.TBATS), std(OTD.STFT.TBATS))
fprintf('|________________|________|________|\n\n')