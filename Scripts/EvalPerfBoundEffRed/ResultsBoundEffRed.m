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
fprintf('|      LSE         |           %.2f        |     %.3e     |\n', CompTime1.LSE , mean(forecastErr1.LSE))
fprintf('|      EDMD        |           %.2f        |     %.3e     |\n', CompTime1.EDMD , mean(forecastErr1.EDMD))
fprintf('|      GPR         |         %.2f        |     %.3e     |\n', CompTime1.GPR , mean(forecastErr1.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('==============Optimal Transport Distance===============\n')
fprintf(' ______________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     |\n')
fprintf('|------------------|-----------|-----------|-----------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.S),mean(OTD1.stft.S),mean(OTD2.rs.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.LSE),mean(OTD1.stft.LSE),mean(OTD2.rs.LSE))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.EDMD),mean(OTD1.stft.EDMD),mean(OTD2.rs.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.GPR),mean(OTD1.stft.GPR),mean(OTD2.rs.GPR))
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
fprintf('|      LSE         |           %.2f        | %.3e | %.3e |\n', CompTime1.LSE , mean(forecastErr1.LSE(~isnan(OTD1.sst.LSE))),median(forecastErr1.LSE(~isnan(OTD1.sst.LSE))) )
fprintf('|      EDMD        |           %.2f        | %.3e | %.3e |\n', CompTime1.EDMD , mean(forecastErr1.EDMD),median(forecastErr1.LSE))
fprintf('|      GPR         |         %.2f        | %.3e | %.3e |\n', CompTime1.GPR , mean(forecastErr1.GPR(~isnan(OTD1.sst.GPR))),median(forecastErr1.GPR(~isnan(OTD1.sst.GPR))))
fprintf('|__________________|_______________________|___________|___________|\n\n')

fprintf('====================Optimal Transport Distance===================\n')
fprintf(' ________________________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     | ConceFT |\n')
fprintf('|------------------|-----------|-----------|-----------|---------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.S),mean(OTD1.stft.S),mean(OTD2.rs.S),mean(OTD3.conceft.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.LSE(~isnan(OTD1.sst.LSE))),mean(OTD1.stft.LSE(~isnan(OTD1.sst.LSE))),mean(OTD2.rs.LSE(~isnan(OTD1.sst.LSE))),mean(OTD3.conceft.LSE(~isnan(OTD3.sst.LSE))))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.EDMD),mean(OTD1.stft.EDMD),mean(OTD2.rs.EDMD),mean(OTD3.conceft.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.GPR(~isnan(OTD1.sst.GPR))),mean(OTD1.stft.GPR(~isnan(OTD1.sst.GPR))),mean(OTD2.rs.GPR(~isnan(OTD1.sst.GPR)))mean(OTD2.conceft.GPR(~isnan(OTD3.sst.GPR))))
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
fprintf('|      LSE         |           %.2f        |     %.3e     |\n', CompTime1.LSE , mean(forecastErr1.LSE))
fprintf('|      EDMD        |           %.2f        |     %.3e     |\n', CompTime1.EDMD , mean(forecastErr1.EDMD))
fprintf('|      GPR         |          %.2f        |     %.3e     |\n', CompTime1.GPR , mean(forecastErr1.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('==============Optimal Transport Distance===============\n')
fprintf(' ______________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     |\n')
fprintf('|------------------|-----------|-----------|-----------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.S),mean(OTD1.stft.S),mean(OTD2.rs.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.LSE),mean(OTD1.stft.LSE),mean(OTD2.rs.LSE))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.EDMD),mean(OTD1.stft.EDMD),mean(OTD2.rs.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.GPR),mean(OTD1.stft.GPR),mean(OTD2.rs.GPR))
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
fprintf('|      LSE         |           %.2f        |     %.3e     |\n', CompTime1.LSE , mean(forecastErr1.LSE))
fprintf('|      EDMD        |          %.2f        |     %.3e     |\n', CompTime1.EDMD , mean(forecastErr1.EDMD))
fprintf('|      GPR         |         %.2f        |     %.3e     |\n', CompTime1.GPR , mean(forecastErr1.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('==============Optimal Transport Distance===============\n')
fprintf(' ______________________________________________________\n')
fprintf('| Extension Method |    SST    |    STFT   |    RS     |\n')
fprintf('|------------------|-----------|-----------|-----------|\n')
fprintf('|  No extension    | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.S),mean(OTD1.stft.S),mean(OTD2.rs.S))
fprintf('|  LSE extension   | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.LSE),mean(OTD1.stft.LSE),mean(OTD2.rs.LSE))
fprintf('|  EDMD extension  | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.EDMD),mean(OTD1.stft.EDMD),mean(OTD2.rs.EDMD))
fprintf('|  GPR extension   | %.3e | %.3e | %.3e |\n', mean(OTD1.sst.GPR),mean(OTD1.stft.GPR),mean(OTD2.rs.GPR))
fprintf('|__________________|___________|___________|___________|\n\n')