%% Display the results of InfluenceSubsigLengthAHM.m and OtherMethAHM.m
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;


%% Display the results of InfluenceSubsigLengthAHM.m

load('../../Results/PerfSubsigLengthAHM') ;

BiasXP = abs(mean(MeanSigExt)) ;
MSEXP = VarSigExt ;
CPUtimeSigExt = mean(CPUtimeXP) ;


fprintf('============= Forecasting Performance =============\n')
fprintf(' ________________________________________________\n')
fprintf('|    Extension    |  Computing  |        MSE      |\n')
fprintf('|     Method      | time (sec.) |   Mean  |   SD  |\n')
fprintf('|-----------------|-------------|---------|-------|\n')

indM = 1 ;
for extM = extMval
    fprintf('| SigExt (M=%4i) |  %8.3f   |  %5.3f  | %5.3f |\n', extM, CPUtimeSigExt(indM) , mean(MSEXP(:,indM)), std(MSEXP(:,indM)) ) ;
    indM = indM + 1 ;
end
indM = indM +1 ;

%% Display the results of OtherMethAHM.m

load('../../Results/PerfOtherMethodsAHM') ;

CPUtimeSYM = mean(CPUtimeXP.SYM) ;
MSESYM = MSE.SYM ;

CPUtimeEDMD = mean(CPUtimeXP.EDMD) ;
MSEEDMD = MSE.EDMD ;

CPUtimeGPR = mean(CPUtimeXP.GPR) ;
MSEGPR = MSE.GPR ;

dataTBATS = table2array( readtable('../../Results/PerfAHM_TBATS.csv','Range','B:D','TreatAsEmpty','NA') ) ;
MSETBATS = dataTBATS(:,2) ;
CPUtimeTBATS = dataTBATS(~isnan(dataTBATS(:,3)),3) ;



fprintf('|  Symmetrization |  %8.3f   |  %5.3f  | %5.3f |\n', CPUtimeSYM , mean(MSESYM), std(MSESYM) ) ;
fprintf('|       EDMD      |  %8.3f   |  %5.3f  | %5.3f |\n', CPUtimeEDMD , mean(MSEEDMD), std(MSEEDMD) ) ;
fprintf('|       GPR       |  %8.3f   |  %5.3f  | %5.3f |\n', CPUtimeGPR , mean(MSEGPR), std(MSEGPR) ) ;
fprintf('|      TBATS      |  %8.3f   |  %5.3f  | %5.3f |\n', CPUtimeTBATS ,mean(MSETBATS), std(MSETBATS))
fprintf('|_________________|_____________|_________|_______|\n\n')