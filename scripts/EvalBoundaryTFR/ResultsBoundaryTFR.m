clear all; close all; clc;


%% PPG
load ../../Results/resultSucForPPGcompMeth ;
fprintf('===============================\n')
fprintf('             PPG               \n')
fprintf('===============================\n')


fprintf('==============SST==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(sstS))
fprintf('|  LSE extension   | %.3e |\n', mean(sstLSE))
fprintf('|__________________|___________|\n\n')

fprintf('=============STFT==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(stftS))
fprintf('|  LSE extension   | %.3e |\n', mean(stftLSE))
fprintf('|__________________|___________|\n\n')

fprintf('==============RS===============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(rsS))
fprintf('|  LSE extension   | %.3e |\n', mean(rsLSE))
fprintf('|__________________|___________|\n\n')

%% THO
load ../../Results/resultSucForTHOcompMeth ;
fprintf('===============================\n')
fprintf('              THO              \n')
fprintf('===============================\n')


fprintf('==============SST==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(sstS))
fprintf('|  LSE extension   | %.3e |\n', mean(sstLSE))
fprintf('|__________________|___________|\n\n')

fprintf('=============STFT==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(stftS))
fprintf('|  LSE extension   | %.3e |\n', mean(stftLSE))
fprintf('|__________________|___________|\n\n')

fprintf('==============RS===============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(rsS))
fprintf('|  LSE extension   | %.3e |\n', mean(rsLSE))
fprintf('|__________________|___________|\n\n')

%% EEG
load ../../Results/resultSucForEEGcompMeth ;
fprintf('===============================\n')
fprintf('              EEG              \n')
fprintf('===============================\n')

fprintf('==============SST==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(sstS))
fprintf('|  LSE extension   | %.3e |\n', mean(sstLSE))
fprintf('|__________________|___________|\n\n')

fprintf('=============STFT==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(stftS))
fprintf('|  LSE extension   | %.3e |\n', mean(stftLSE))
fprintf('|__________________|___________|\n\n')

fprintf('===============RS==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(rsS))
fprintf('|  LSE extension   | %.3e |\n', mean(rsLSE))
fprintf('|__________________|___________|\n\n')

%% ECG
load ../../Results/resultSucForECGcompMeth ;
fprintf('===============================\n')
fprintf('              ECG              \n')
fprintf('===============================\n')

fprintf('==============SST==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(sstS))
fprintf('|  LSE extension   | %.3e |\n', mean(sstLSE))
fprintf('|__________________|___________|\n\n')

fprintf('=============STFT==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(stftS))
fprintf('|  LSE extension   | %.3e |\n', mean(stftLSE))
fprintf('|__________________|___________|\n\n')

fprintf('===============RS==============\n')
fprintf(' ______________________________\n')
fprintf('| Extension Method |    OTD    | \n')
fprintf('|------------------|-----------|\n')
fprintf('|  Short signal    | %.3e |\n', mean(rsS))
fprintf('|  LSE extension   | %.3e |\n', mean(rsLSE))
fprintf('|__________________|___________|\n\n')