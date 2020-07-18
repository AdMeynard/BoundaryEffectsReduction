function [MultiTaperSTFT, tfrtic, MultiTaper, MultiTaperAll, tfrrstic] = rsSTFTmultitaper(x, lowFreq, highFreq, alpha, WinLen, dim, supp)

%
% Usage: 
% 	[tfrsq, ConceFT, tfrsqtic] = sqSTFTmultitaper(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp)
%
% MT = 1: ordinary SST; MT > 1: ConceFT
% alpha: resolution in the frequency axis
% WinLen, dim, supp: for hermf.m
%
% Example:
% 	[tfrsq, MultiTaper, tfrsqtic] = sqSTFT([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6);


N = length(x) ;

%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       	%% generate the window for short time Fourier transform (STFT)
[h, Dh, ~] = hermf(WinLen, dim, supp) ;


%=======================================

[tfr, tfrtic , tfrrs, tfrrstic] = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, h(1,:)', Dh(1,:)', 1);

MultiTaperAll = zeros(size(tfrrs,1), size(tfrrs,2), dim) ;
MultiTaperAll(:,:,1) = tfrrs ;
MultiTaperSTFT = tfr ;

fprintf(['MultiTaper total: ',num2str(dim),'; now:     ']) ;
for ii = 2: dim
	fprintf('\b\b\b\b') ;	tmp = sprintf('%4d',ii) ; fprintf([tmp]) ;
	[tfr, ~, tfrrs, tfrrstic] = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, h(ii,:)', Dh(ii,:)', 1);

 	MultiTaperAll(:, :, ii) = tfrrs ;
	MultiTaperSTFT = MultiTaperSTFT + tfr ;
end

MultiTaper = mean(MultiTaperAll, 3) ;
MultiTaperSTFT = MultiTaperSTFT / dim ;
fprintf('\n') ;

end
