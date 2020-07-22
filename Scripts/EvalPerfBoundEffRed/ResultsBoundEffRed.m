%% Performance of BoundEffRed on PPG,THO, ECG, and EEG (Tables II and II in the paper)
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;

%% PPG
load ../../Results/resultSucForPPGcompMeth ;
fprintf('================================================================\n')
fprintf('                              PPG                               \n')
fprintf('================================================================\n')

fprintf('====================Forecasting Performance=====================\n')
fprintf(' _______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.)  | Forecasting Error |\n')
fprintf('|------------------|------------------------|-------------------|\n')
fprintf('|     SigExt       |           %.3f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE) )
fprintf('|  Symmetrization  |           %.3f        |     %.3e     |\n', CompTime.SYM , mean(forecastErr.SYM) )
fprintf('|      EDMD        |           %.3f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD) )
fprintf('|      GPR         |         %.3f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR) )
fprintf('|__________________|________________________|___________________|\n\n')

fprintf('=========Optimal Transport Distance=========\n')
fprintf(' __________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |ConceFT|\n')
fprintf('|------------------|-------|-------|-------|\n')
fprintf('|      SigExt      | %.3f | %.3f | %.3f |\n', mean(OTD.sst.LSE./OTD.sst.S),mean(OTD.stft.LSE./OTD.stft.S),mean(OTD.conceft.LSE./OTD.conceft.S))
fprintf('|  Symmetrization  | %.3f | %.3f | %.3f |\n', mean(OTD.sst.SYM./OTD.sst.S),mean(OTD.stft.SYM./OTD.stft.S),mean(OTD.conceft.SYM./OTD.conceft.S))
fprintf('|  EDMD extension  | %.3f | %.3f | %.3f |\n', mean(OTD.sst.EDMD./OTD.sst.S),mean(OTD.stft.EDMD./OTD.stft.S),mean(OTD.conceft.EDMD./OTD.conceft.S))
fprintf('|  GPR extension   | %.3f | %.3f | %.3f |\n', mean(OTD.sst.GPR./OTD.sst.S),mean(OTD.stft.GPR./OTD.stft.S),mean(OTD.conceft.GPR./OTD.conceft.S))
fprintf('|__________________|_______|_______|_______|\n\n')

%% THO
load ../../Results/resultSucForTHOcompMeth ;
fprintf('====================================================================\n')
fprintf('                              THO                                  \n')
fprintf('====================================================================\n')

fprintf('============================Forecasting Performance==================\n')
fprintf(' ___________________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.)  |   Forecasting Error   |\n')
fprintf('|                  |                        |    Mean   |   Median  |\n')
fprintf('|------------------|------------------------|-----------|-----------|\n')
fprintf('|     SigExt       |           %.3f        | %.3e | %.3e |\n', CompTime.LSE , mean(forecastErr.LSE(~isnan(OTD.sst.LSE))),median(forecastErr.LSE(~isnan(OTD.sst.LSE))) )
fprintf('|  Symmetrization  |           %.3f        | %.3e | %.3e |\n', CompTime.SYM , mean(forecastErr.SYM), median(forecastErr.SYM) )
fprintf('|      EDMD        |          %.3f        | %.3e | %.3e |\n', CompTime.EDMD , mean(forecastErr.EDMD), median(forecastErr.EDMD) )
fprintf('|      GPR         |         %.3f        | %.3e | %.3e |\n', CompTime.GPR , mean(forecastErr.GPR(~isnan(OTD.sst.GPR))), median(forecastErr.GPR(~isnan(OTD.sst.GPR))) )
fprintf('|__________________|________________________|___________|___________|\n\n')

figure;
subplot(2,1,1) ;
histogram( log10 ( forecastErr.LSE(~isnan(OTD.sst.LSE)) ) , 100)

meanH = mean(forecastErr.LSE(~isnan(OTD.sst.LSE))) ;
medH = median(forecastErr.LSE(~isnan(OTD.sst.LSE))) ;
xline(log10(meanH),'r-',{'Mean'},'linewidth',2,'fontsize',16);
xline(log10(medH),'r-',{'Median'},'linewidth',2,'fontsize',16);

u = xticks ; str = compose('10^{%.1f}',u) ;
xticklabels(str); grid on;
xlabel('MSE'); ylabel('Histogram');
set(gca,'fontsize',20)

load ../../Results/FailureCase ;
N = length(x);
t = linspace(0,(N-1)/fs,N);
tt = linspace(-L/fs,(N-1)/fs + L/fs,N+2*L);

subplot(2,1,2);
plot(tt,xxLSE,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
legend({'{\sf SigExt} Extended signal','Ground truth Extended signal','Original signal'},'interpreter','latex'); 
xlabel('Time (s)'); ylabel('Signals'); xlim([40 tt(end)]); ylim([-20 20]);
set(gca,'fontsize',20)


fprintf('============Optimal Transport Distance=============\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |  RS   |ConceFT|\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %.3f | %.3f | %.3f | %.3f |\n', mean(OTD.sst.LSE(~isnan(OTD.sst.LSE))./OTD.sst.S(~isnan(OTD.sst.LSE))),mean(OTD.stft.LSE(~isnan(OTD.sst.LSE))./OTD.stft.S(~isnan(OTD.sst.LSE))),mean(OTD.rs.LSE(~isnan(OTD.sst.LSE))./OTD.rs.S(~isnan(OTD.sst.LSE))),mean(OTD.conceft.LSE(~isnan(OTD.sst.LSE))./OTD.conceft.S(~isnan(OTD.sst.LSE))))
fprintf('|  Symmetrization  | %.3f | %.3f | %.3f | %.3f |\n', mean(OTD.sst.SYM./OTD.sst.S),mean(OTD.stft.SYM./OTD.stft.S),mean(OTD.rs.SYM./OTD.rs.S),mean(OTD.conceft.SYM./OTD.conceft.S))
fprintf('|  EDMD extension  | %.3f | %.3f | %.3f | %.3f |\n', mean(OTD.sst.EDMD./OTD.sst.S),mean(OTD.stft.EDMD./OTD.stft.S),mean(OTD.rs.EDMD./OTD.rs.S),mean(OTD.conceft.EDMD./OTD.conceft.S))
fprintf('|  GPR extension   | %.3f | %.3f | %.3f | %.3f |\n', mean(OTD.sst.GPR(~isnan(OTD.sst.GPR))./OTD.sst.S(~isnan(OTD.sst.GPR))),mean(OTD.stft.GPR(~isnan(OTD.sst.GPR))./OTD.stft.S(~isnan(OTD.sst.GPR))),mean(OTD.rs.GPR(~isnan(OTD.sst.GPR))./OTD.rs.S(~isnan(OTD.sst.GPR))),mean(OTD.conceft.GPR(~isnan(OTD.sst.GPR))./OTD.conceft.S(~isnan(OTD.sst.GPR))) )
fprintf('|__________________|_______|_______|_______|_______|\n\n')

%% EEG
load ../../Results/resultSucForEEGcompMeth ;
fprintf('================================================================\n')
fprintf('                             EEG                                \n')
fprintf('================================================================\n')

fprintf('===================Forecasting Performance======================\n')
fprintf(' ______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|------------------|-----------------------|-------------------|\n')
fprintf('|     SigExt       |           %.2f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE))
fprintf('|      EDMD        |           %.2f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD))
fprintf('|      GPR         |          %.2f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('=========Optimal Transport Distance========\n')
fprintf(' __________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |  RS   |\n')
fprintf('|------------------|-------|-------|-------|\n')
fprintf('|      SigExt      | %.3f | %.3f | %.3f |\n', mean(OTD.sst.LSE./OTD.sst.S),mean(OTD.stft.LSE./OTD.stft.S),mean(OTD.rs.LSE./OTD.rs.S))
%fprintf('|  Symmetrization  | %.3f | %.3f | %.3f |\n', mean(OTD.sst.SYM./OTD.sst.S),mean(OTD.stft.SYM./OTD.stft.S),mean(OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  | %.3f | %.3f | %.3f |\n', mean(OTD.sst.EDMD./OTD.sst.S),mean(OTD.stft.EDMD./OTD.stft.S),mean(OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   | %.3f | %.3f | %.3f |\n', mean(OTD.sst.GPR./OTD.sst.S),mean(OTD.stft.GPR./OTD.stft.S),mean(OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_______|_______|_______|\n\n')

%% ECG
load ../../Results/resultSucForECGcompMeth ;
fprintf('================================================================\n')
fprintf('                            ECG                                 \n')
fprintf('================================================================\n')

fprintf('===================Forecasting Performance======================\n')
fprintf(' ______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|------------------|-----------------------|-------------------|\n')
fprintf('|     SigExt       |           %.2f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE))
fprintf('|      EDMD        |          %.2f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD))
fprintf('|      GPR         |         %.2f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('=========Optimal Transport Distance========\n')
fprintf(' __________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |  RS   |\n')
fprintf('|------------------|-------|-------|-------|\n')
fprintf('|     SigExt       | %.3f | %.3f | %.3f |\n', mean(OTD.sst.LSE./OTD.sst.S),mean(OTD.stft.LSE./OTD.stft.S),mean(OTD.rs.LSE./OTD.rs.S))
%fprintf('|  Symmetrization  | %.3f | %.3f | %.3f |\n', mean(OTD.sst.SYM./OTD.sst.S),mean(OTD.stft.SYM./OTD.stft.S),mean(OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  | %.3f | %.3f | %.3f |\n', mean(OTD.sst.EDMD./OTD.sst.S),mean(OTD.stft.EDMD./OTD.stft.S),mean(OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   | %.3f | %.3f | %.3f |\n', mean(OTD.sst.GPR./OTD.sst.S),mean(OTD.stft.GPR./OTD.stft.S),mean(OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_______|_______|_______|\n\n')