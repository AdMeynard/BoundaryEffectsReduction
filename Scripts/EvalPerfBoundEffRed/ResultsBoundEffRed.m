%% Performance of BoundEffRed on PPG,THO, ECG, and EEG (Tables II and II in the paper)
% Author: Adrien MEYNARD
% Email: adrien.meynard@duke.edu

clear all; close all; clc;

%% PPG
load ../../Results/resultSucForPPGcompMeth ;
fprintf('================================================================\n')
fprintf('                              PPG                               \n')
fprintf('================================================================\n')

fprintf('===================== Forecasting Performance ==================\n')
fprintf(' _______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|                  |                       |   Mean  |    SD   |\n')
fprintf('|------------------|-----------------------|---------|---------|\n')
fprintf('|     SigExt       |         %7.3f       |  %.3f  |  %.3f  |\n', CompTime.LSE , mean(forecastErr.LSE(~isnan(OTD.sst.LSE))),std(forecastErr.LSE(~isnan(OTD.sst.LSE))) )
fprintf('|  Symmetrization  |         %7.3f       |  %.3f  |  %.3f  |\n', CompTime.SYM , mean(forecastErr.SYM), std(forecastErr.SYM) )
fprintf('|      EDMD        |         %7.3f       |  %.3f  |  %.3f  |\n', CompTime.EDMD , mean(forecastErr.EDMD), std(forecastErr.EDMD) )
fprintf('|      GPR         |         %7.3f       |  %.3f  |  %.3f  |\n', CompTime.GPR , mean(forecastErr.GPR(~isnan(OTD.sst.GPR))), std(forecastErr.GPR(~isnan(OTD.sst.GPR))) )
fprintf('|__________________|_______________________|_________|_________|\n\n')

fprintf('========== Averaged Performance Index D ===========\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |ConceFT|   RS  |\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.LSE./OTD.sst.S),mean(OTD.stft.LSE./OTD.stft.S),mean(OTD.conceft.LSE./OTD.conceft.S),mean(OTD.rs.LSE./OTD.rs.S))
fprintf('|  Symmetrization  | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.SYM./OTD.sst.S),mean(OTD.stft.SYM./OTD.stft.S),mean(OTD.conceft.SYM./OTD.conceft.S),mean(OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.EDMD./OTD.sst.S),mean(OTD.stft.EDMD./OTD.stft.S),mean(OTD.conceft.EDMD./OTD.conceft.S),mean(OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.GPR./OTD.sst.S),mean(OTD.stft.GPR./OTD.stft.S),mean(OTD.conceft.GPR./OTD.conceft.S),mean(OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_______|_______|_______|_______|\n\n')

fprintf('========== SD of the Performance Index D ===========\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |  RS   |ConceFT|\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.LSE./OTD.sst.S),std(OTD.stft.LSE./OTD.stft.S),std(OTD.rs.LSE./OTD.rs.S),std(OTD.conceft.LSE./OTD.conceft.S))
fprintf('|  Symmetrization  | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.SYM./OTD.sst.S),std(OTD.stft.SYM./OTD.stft.S),std(OTD.rs.SYM./OTD.rs.S),std(OTD.conceft.SYM./OTD.conceft.S))
fprintf('|  EDMD extension  | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.EDMD./OTD.sst.S),std(OTD.stft.EDMD./OTD.stft.S),std(OTD.rs.EDMD./OTD.rs.S),std(OTD.conceft.EDMD./OTD.conceft.S))
fprintf('|  GPR extension   | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.GPR./OTD.sst.S),std(OTD.stft.GPR./OTD.stft.S),std(OTD.rs.GPR./OTD.rs.S),std(OTD.conceft.GPR./OTD.conceft.S) )
fprintf('|__________________|_______|_______|_______|_______|\n\n')

fprintf('============= T-Test vs. SigExt ===========\n')
fprintf(' __________________________________________\n')
fprintf('| Extension Method | SST |STFT|ConceFT| RS |\n')
fprintf('|------------------|-----|----|-------|----|\n')
fprintf('|  Symmetrization  |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.SYM./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.SYM./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.SYM./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.EDMD./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.EDMD./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.EDMD./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.GPR./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.GPR./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.GPR./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_____|____|_______|____|\n\n')

%% THO
load ../../Results/resultSucForTHOcompMeth ;
fprintf('====================================================================\n')
fprintf('                              THO                                  \n')
fprintf('====================================================================\n')

fprintf('======================= Forecasting Performance ==================\n')
fprintf(' _______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|                  |                       |   Mean  |    SD   |\n')
fprintf('|------------------|-----------------------|---------|---------|\n')
fprintf('|     SigExt       |         %7.3f       |  %5.3f  |  %5.3f  |\n', CompTime.LSE , mean(forecastErr.LSE(~isnan(OTD.sst.LSE))),std(forecastErr.LSE(~isnan(OTD.sst.LSE))) )
fprintf('|  Symmetrization  |         %7.3f       |  %5.3f  |  %5.3f  |\n', CompTime.SYM , mean(forecastErr.SYM), std(forecastErr.SYM) )
fprintf('|      EDMD        |         %7.3f       |  %5.3f  |  %5.3f  |\n', CompTime.EDMD , mean(forecastErr.EDMD), std(forecastErr.EDMD) )
fprintf('|      GPR         |         %7.3f       |  %5.3f  |  %5.3f  |\n', CompTime.GPR , mean(forecastErr.GPR(~isnan(OTD.sst.GPR))), std(forecastErr.GPR(~isnan(OTD.sst.GPR))) )
fprintf('|__________________|_______________________|_________|_________|\n\n')

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


fprintf('========== Averaged Performance Index D ===========\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |  RS   |ConceFT|\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.LSE(~isnan(OTD.sst.LSE))./OTD.sst.S(~isnan(OTD.sst.LSE))),mean(OTD.stft.LSE(~isnan(OTD.sst.LSE))./OTD.stft.S(~isnan(OTD.sst.LSE))),mean(OTD.rs.LSE(~isnan(OTD.sst.LSE))./OTD.rs.S(~isnan(OTD.sst.LSE))),mean(OTD.conceft.LSE(~isnan(OTD.sst.LSE))./OTD.conceft.S(~isnan(OTD.sst.LSE))))
fprintf('|  Symmetrization  | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.SYM./OTD.sst.S),mean(OTD.stft.SYM./OTD.stft.S),mean(OTD.rs.SYM./OTD.rs.S),mean(OTD.conceft.SYM./OTD.conceft.S))
fprintf('|  EDMD extension  | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.EDMD./OTD.sst.S),mean(OTD.stft.EDMD./OTD.stft.S),mean(OTD.rs.EDMD./OTD.rs.S),mean(OTD.conceft.EDMD./OTD.conceft.S))
fprintf('|  GPR extension   | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.GPR(~isnan(OTD.sst.GPR))./OTD.sst.S(~isnan(OTD.sst.GPR))),mean(OTD.stft.GPR(~isnan(OTD.sst.GPR))./OTD.stft.S(~isnan(OTD.sst.GPR))),mean(OTD.rs.GPR(~isnan(OTD.sst.GPR))./OTD.rs.S(~isnan(OTD.sst.GPR))),mean(OTD.conceft.GPR(~isnan(OTD.sst.GPR))./OTD.conceft.S(~isnan(OTD.sst.GPR))) )
fprintf('|__________________|_______|_______|_______|_______|\n\n')

fprintf('========== SD of the Performance Index D ===========\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |  RS   |ConceFT|\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.LSE(~isnan(OTD.sst.LSE))./OTD.sst.S(~isnan(OTD.sst.LSE))),std(OTD.stft.LSE(~isnan(OTD.sst.LSE))./OTD.stft.S(~isnan(OTD.sst.LSE))),std(OTD.rs.LSE(~isnan(OTD.sst.LSE))./OTD.rs.S(~isnan(OTD.sst.LSE))),std(OTD.conceft.LSE(~isnan(OTD.sst.LSE))./OTD.conceft.S(~isnan(OTD.sst.LSE))))
fprintf('|  Symmetrization  | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.SYM./OTD.sst.S),std(OTD.stft.SYM./OTD.stft.S),std(OTD.rs.SYM./OTD.rs.S),std(OTD.conceft.SYM./OTD.conceft.S))
fprintf('|  EDMD extension  | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.EDMD./OTD.sst.S),std(OTD.stft.EDMD./OTD.stft.S),std(OTD.rs.EDMD./OTD.rs.S),std(OTD.conceft.EDMD./OTD.conceft.S))
fprintf('|  GPR extension   | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.GPR(~isnan(OTD.sst.GPR))./OTD.sst.S(~isnan(OTD.sst.GPR))),std(OTD.stft.GPR(~isnan(OTD.sst.GPR))./OTD.stft.S(~isnan(OTD.sst.GPR))),std(OTD.rs.GPR(~isnan(OTD.sst.GPR))./OTD.rs.S(~isnan(OTD.sst.GPR))),std(OTD.conceft.GPR(~isnan(OTD.sst.GPR))./OTD.conceft.S(~isnan(OTD.sst.GPR))) )
fprintf('|__________________|_______|_______|_______|_______|\n\n')

fprintf('============= T-Test vs. SigExt ===========\n')
fprintf(' __________________________________________\n')
fprintf('| Extension Method | SST |STFT|ConceFT| RS |\n')
fprintf('|------------------|-----|----|-------|----|\n')
fprintf('|  Symmetrization  |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.SYM./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.SYM./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.SYM./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.EDMD./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.EDMD./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.EDMD./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.GPR./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.GPR./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.GPR./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_____|____|_______|____|\n\n')

%% EEG
load ../../Results/resultSucForEEGcompMeth ;
fprintf('================================================================\n')
fprintf('                             EEG                                \n')
fprintf('================================================================\n')

fprintf('=================== Forecasting Performance ===================\n')
fprintf(' ______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|------------------|-----------------------|-------------------|\n')
fprintf('|     SigExt       |          %5.2f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE))
fprintf('|      EDMD        |          %5.2f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD))
fprintf('|      GPR         |          %5.2f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('========== Median Performance Index D =============\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |ConceFT|   RS  |\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %5.3f | %5.3f | %5.3f | %5.3f |\n', median(OTD.sst.LSE./OTD.sst.S),median(OTD.stft.LSE./OTD.stft.S),median(OTD.conceft.LSE./OTD.conceft.S),median(OTD.rs.LSE./OTD.rs.S))
fprintf('|  Symmetrization  | %5.3f | %5.3f | %5.3f | %5.3f |\n', median(OTD.sst.SYM./OTD.sst.S),median(OTD.stft.SYM./OTD.stft.S),median(OTD.conceft.SYM./OTD.conceft.S),median(OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  | %5.3f | %5.3f | %5.3f | %5.3f |\n', median(OTD.sst.EDMD./OTD.sst.S),median(OTD.stft.EDMD./OTD.stft.S),median(OTD.conceft.EDMD./OTD.conceft.S),median(OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   | %5.3f | %5.3f | %5.3f | %5.3f |\n', median(OTD.sst.GPR./OTD.sst.S),median(OTD.stft.GPR./OTD.stft.S),median(OTD.conceft.GPR./OTD.conceft.S),median(OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_______|_______|_______|_______|\n\n')

%% ECG
load ../../Results/resultSucForECGcompMeth ;
fprintf('================================================================\n')
fprintf('                            ECG                                 \n')
fprintf('================================================================\n')

fprintf('=================== Forecasting Performance ===================\n')
fprintf(' ______________________________________________________________\n')
fprintf('| Extension Method | Computing time (sec.) | Forecasting Error |\n')
fprintf('|------------------|-----------------------|-------------------|\n')
fprintf('|     SigExt       |         %6.2f        |     %.3e     |\n', CompTime.LSE , mean(forecastErr.LSE))
fprintf('|      EDMD        |         %6.2f        |     %.3e     |\n', CompTime.EDMD , mean(forecastErr.EDMD))
fprintf('|      GPR         |         %6.2f        |     %.3e     |\n', CompTime.GPR , mean(forecastErr.GPR))
fprintf('|__________________|_______________________|___________________|\n\n')

fprintf('========== Averaged Performance Index D ===========\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |ConceFT|   RS  |\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.LSE./OTD.sst.S),mean(OTD.stft.LSE./OTD.stft.S),mean(OTD.conceft.LSE./OTD.conceft.S),mean(OTD.rs.LSE./OTD.rs.S))
fprintf('|  Symmetrization  | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.SYM./OTD.sst.S),mean(OTD.stft.SYM./OTD.stft.S),mean(OTD.conceft.SYM./OTD.conceft.S),mean(OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.EDMD./OTD.sst.S),mean(OTD.stft.EDMD./OTD.stft.S),mean(OTD.conceft.EDMD./OTD.conceft.S),mean(OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   | %5.3f | %5.3f | %5.3f | %5.3f |\n', mean(OTD.sst.GPR./OTD.sst.S),mean(OTD.stft.GPR./OTD.stft.S),mean(OTD.conceft.GPR./OTD.conceft.S),mean(OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_______|_______|_______|_______|\n\n')

fprintf('========== SD of the Performance Index D ===========\n')
fprintf(' __________________________________________________\n')
fprintf('| Extension Method |  SST  |  STFT |  RS   |ConceFT|\n')
fprintf('|------------------|-------|-------|-------|-------|\n')
fprintf('|      SigExt      | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.LSE./OTD.sst.S),std(OTD.stft.LSE./OTD.stft.S),std(OTD.rs.LSE./OTD.rs.S),std(OTD.conceft.LSE./OTD.conceft.S))
fprintf('|  Symmetrization  | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.SYM./OTD.sst.S),std(OTD.stft.SYM./OTD.stft.S),std(OTD.rs.SYM./OTD.rs.S),std(OTD.conceft.SYM./OTD.conceft.S))
fprintf('|  EDMD extension  | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.EDMD./OTD.sst.S),std(OTD.stft.EDMD./OTD.stft.S),std(OTD.rs.EDMD./OTD.rs.S),std(OTD.conceft.EDMD./OTD.conceft.S))
fprintf('|  GPR extension   | %5.3f | %5.3f | %5.3f | %5.3f |\n', std(OTD.sst.GPR./OTD.sst.S),std(OTD.stft.GPR./OTD.stft.S),std(OTD.rs.GPR./OTD.rs.S),std(OTD.conceft.GPR./OTD.conceft.S) )
fprintf('|__________________|_______|_______|_______|_______|\n\n')

fprintf('============= T-Test vs. SigExt ===========\n')
fprintf(' __________________________________________\n')
fprintf('| Extension Method | SST |STFT|ConceFT| RS |\n')
fprintf('|------------------|-----|----|-------|----|\n')
fprintf('|  Symmetrization  |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.SYM./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.SYM./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.SYM./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.SYM./OTD.rs.S))
fprintf('|  EDMD extension  |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.EDMD./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.EDMD./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.EDMD./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.EDMD./OTD.rs.S))
fprintf('|  GPR extension   |  %u  |  %u |   %u   |  %u |\n', ttest(OTD.sst.LSE./OTD.sst.S,OTD.sst.GPR./OTD.sst.S),ttest(OTD.stft.LSE./OTD.stft.S,OTD.stft.GPR./OTD.stft.S),ttest(OTD.conceft.LSE./OTD.conceft.S,OTD.conceft.GPR./OTD.conceft.S),ttest(OTD.rs.LSE./OTD.rs.S,OTD.rs.GPR./OTD.rs.S))
fprintf('|__________________|_____|____|_______|____|\n\n')