clear all; close all; clc;

%% PPG
load ../../Results/resultSucForPPGcompMeth ;
fprintf('================================================================\n')
fprintf('                              PPG                               \n')
fprintf('================================================================\n')

fprintf('===================Forecasting Performance======================\n')
fprintf(' ______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|------------------|-----------------------|-------------------|\n')
fprintf('|      LSE         |           %.2f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE))
fprintf('|      EDMD        |           %.2f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD))
fprintf('|      GPR         |         %.2f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('==============Optimal Transport Distance===============\n')
fprintf(' ______________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     |\n')
fprintf('|------------------|-----------|-----------|-----------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e |\n', mean(OTD.sst.S),mean(OTD.stft.S),mean(OTD.rs.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e |\n', mean(OTD.sst.LSE),mean(OTD.stft.LSE),mean(OTD.rs.LSE))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e |\n', mean(OTD.sst.EDMD),mean(OTD.stft.EDMD),mean(OTD.rs.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e |\n', mean(OTD.sst.GPR),mean(OTD.stft.GPR),mean(OTD.rs.GPR))
fprintf('|__________________|___________|___________|___________|\n\n')

%% THO
load ../../Results/resultSucForTHOcompMeth ;
fprintf('====================================================================\n')
fprintf('                              THO                                  \n')
fprintf('====================================================================\n')

fprintf('===========================Forecasting Performance==================\n')
fprintf(' __________________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) |   Forecasting Error   |\n')
fprintf('|                  |                       |    Mean   |   Median  |\n')
fprintf('|------------------|-----------------------|-----------|-----------|\n')
fprintf('|      LSE         |           %.2f        | %.3e | %.3e |\n', CompTime.LSE , mean(forecastErr.LSE(~isnan(OTD.sst.LSE))),median(forecastErr.LSE(~isnan(OTD.sst.LSE))) )
fprintf('|      EDMD        |           %.2f        | %.3e | %.3e |\n', CompTime.EDMD , mean(forecastErr.EDMD),median(forecastErr.LSE))
fprintf('|      GPR         |         %.2f        | %.3e | %.3e |\n', CompTime.GPR , mean(forecastErr.GPR(~isnan(OTD.sst.GPR))),median(forecastErr.GPR(~isnan(OTD.sst.GPR))))
fprintf('|__________________|_______________________|___________|___________|\n\n')

fprintf('====================Optimal Transport Distance===================\n')
fprintf(' ________________________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     | ConceFT |\n')
fprintf('|------------------|-----------|-----------|-----------|---------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD.sst.S),mean(OTD.stft.S),mean(OTD.rs.S),mean(OTD.conceft.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD.sst.LSE(~isnan(OTD.sst.LSE))),mean(OTD.stft.LSE(~isnan(OTD.sst.LSE))),mean(OTD.rs.LSE(~isnan(OTD.sst.LSE))),mean(OTD.conceft.LSE(~isnan(OTD.sst.LSE))))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD.sst.EDMD),mean(OTD.stft.EDMD),mean(OTD.rs.EDMD),mean(OTD.conceft.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD.sst.GPR(~isnan(OTD.sst.GPR))),mean(OTD.stft.GPR(~isnan(OTD.sst.GPR))),mean(OTD.rs.GPR(~isnan(OTD.sst.GPR))),mean(OTD.conceft.GPR(~isnan(OTD.sst.GPR))))
fprintf('|__________________|___________|___________|___________|_________|\n\n')

%% EEG
load ../../Results/resultSucForEEGcompMeth ;
fprintf('================================================================\n')
fprintf('                             EEG                                \n')
fprintf('================================================================\n')

fprintf('===================Forecasting Performance======================\n')
fprintf(' ______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|------------------|-----------------------|-------------------|\n')
fprintf('|      LSE         |           %.2f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE))
fprintf('|      EDMD        |           %.2f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD))
fprintf('|      GPR         |          %.2f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('==============Optimal Transport Distance===============\n')
fprintf(' ______________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     |\n')
fprintf('|------------------|-----------|-----------|-----------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e |\n', mean(OTD.sst.S),mean(OTD.stft.S),mean(OTD.rs.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e |\n', mean(OTD.sst.LSE),mean(OTD.stft.LSE),mean(OTD.rs.LSE))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e |\n', mean(OTD.sst.EDMD),mean(OTD.stft.EDMD),mean(OTD.rs.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e |\n', mean(OTD.sst.GPR),mean(OTD.stft.GPR),mean(OTD.rs.GPR))
fprintf('|__________________|___________|___________|___________|\n\n')

%% ECG
load ../../Results/resultSucForECGcompMeth ;
fprintf('================================================================\n')
fprintf('                            ECG                                 \n')
fprintf('================================================================\n')

fprintf('===================Forecasting Performance======================\n')
fprintf(' ______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|------------------|-----------------------|-------------------|\n')
fprintf('|      LSE         |           %.2f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE))
fprintf('|      EDMD        |          %.2f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD))
fprintf('|      GPR         |         %.2f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('==============Optimal Transport Distance===============\n')
fprintf(' ______________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     |\n')
fprintf('|------------------|-----------|-----------|-----------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e |\n', mean(OTD.sst.S),mean(OTD.stft.S),mean(OTD.rs.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e |\n', mean(OTD.sst.LSE),mean(OTD.stft.LSE),mean(OTD.rs.LSE))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e |\n', mean(OTD.sst.EDMD),mean(OTD.stft.EDMD),mean(OTD.rs.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e |\n', mean(OTD.sst.GPR),mean(OTD.stft.GPR),mean(OTD.rs.GPR))
fprintf('|__________________|___________|___________|___________|\n\n')