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
    fprintf('| SigExt (M=%4i) |     %.3f   |  %.3f  | %.3f |\n', extM, CPUtimeSigExt(indM) , mean(MSEXP(:,indM)), std(MSEXP(:,indM)) ) ;
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


fprintf('|  Symmetrization |     %.3f   |  %.3f  | %.3f |\n', CPUtimeSYM , mean(MSESYM), std(MSESYM) ) ;
fprintf('|       EDMD      |     %.3f   |  %.3f  | %.3f |\n', CPUtimeEDMD , mean(MSEEDMD), std(MSEEDMD) ) ;
fprintf('|       GPR       |   %.3f   |  %.3f  | %.3f |\n', CPUtimeGPR , mean(MSEGPR), std(MSEGPR) ) ;

fprintf('|_________________|_____________|_________|_______|\n\n')