function [tfr, tfrsq, ConceFT, tfrtic] = ConceFT_sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Second, Smooth)
% CONCEFT_SQSTFT_C ConceFT Synchrosqueezing of the STFT
% Usage: 
% 	[tfr, tfrsq, ConceFT, tfrtic] = ConceFT_sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Second, Smooth)
%
% Input:
%   x: signal to be analized
%   lowFreq: relative lower frequency (0<lowFreq<Highfreq)
%   highFreq: relative higher frequency (lowFreq<Highfreq<0.5)
%   alpha: resolution in the frequency axis
%   hop: hop size (in samples)
%   WinLen: Window length (in samples)
%   dim: order of the Hermite function used as window
%   supp: half time support of the Hermite function used as window
%   MT = 1: ordinary RS; MT > 1: ConceFT
%   Second: if 1, compute the second-order synchrosqueezing
%   Smooth: if 1, compute the smoothed version of the synchrosqueezing
% 
% Output:
%   tfr: STFT
%   tfrsq: Synchrosqueezing of STFT
%   ConceFT: ConceFT synchrosqueezing of STFT
%   tfrtic: frequencies for which TFRs are evaluated

method = 'MEAN' ;

%% Ordinary SST
[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

if ~Second
    [tfr, tfrsq, tfrtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', Smooth);
else
    [tfr, ~, tfrsq, tfrtic] = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', dwindow(Dh(1,:)'));
end

%% Multitapering

ConceFT = abs(tfrsq) ;
if strcmp(method, 'MEAN')
else
	ConceFTcum = zeros([size(tfrsq) MT+1]) ;
	ConceFTcum(:,:,1) = tfrsq ;
end

if MT > 1

		%% Conceft

    for ii = 1: MT
		rv = randn(1, dim) + sqrt(-1)*randn(1, dim) ; rv = rv ./ norm(rv) ;
		rh = rv * h ; 
		rDh = rv * Dh ;

		if ~Second
			[~, tfrsqX, tfrtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, rh', rDh', Smooth);
		else
			[~, ~, tfrsqX, tfrtic] = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, rh', rDh', dwindow(rDh'));
		end

		if strcmp(method, 'MEAN')
	 		ConceFT = ConceFT + abs(tfrsqX) ;
		else
			ConceFTcum(:,:,ii+1) = tfrsqX ;
		end
    end

	if strcmp(method, 'MEAN')
    	ConceFT = ConceFT ./ (MT+1) ;
    else
        for aa = 1:size(ConceFTcum,1)
            for bb = 1:size(ConceFTcum,2)
                tmp = squeeze(ConceFTcum(aa,bb,:)) ; 
                tmp2 = quantile(abs(tmp), .98) ;
                tmp(find(abs(tmp)>tmp2)) = 0 ; %tmp2.*tmp(find(abs(tmp)>tmp2))./abs(tmp(find(abs(tmp)>tmp2))) ;
                ConceFTcum(aa, bb, :) = tmp ;
            end
        end
        ConceFT = median(ConceFTcum, 3) ;
	end

end
