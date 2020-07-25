function [tfr, tfrrs, ConceFT, tfrtic] = ConceFT_rsSTFT_RT(x, lowFreq, highFreq, alpha, hop, WinLen, LL, dim, supp, MT)
% CONCEFT_RSSTFT ConceFT Reassigment of the STFT
% Usage: 
% 	[tfr, tfrrs, ConceFT, tfrtic] = ConceFT_rsSTFT(x, lowFreq, highFreq, alpha, hop, WinLen, LL, dim, supp, MT)
%
% Input:
%   x: signal to be analized
%   lowFreq: relative lower frequency (0<lowFreq<Highfreq)
%   highFreq: relative higher frequency (lowFreq<Highfreq<0.5)
%   alpha: resolution in the frequency axis
%   hop: hop size (in samples)
%   WinLen: Window length (in samples)
%   LL: Wextension length
%   dim: order of the Hermite function used as window
%   supp: half time support of the Hermite function used as window 
%   MT = 1: ordinary RS; MT > 1: ConceFT
% 
% Output:
%   tfr: STFT
%   tfrrs: Reassigment of STFT
%   ConceFT: ConceFT RS of STFT
%   tfrtic: frequencies for which TFRs are evaluated


%% Ordinary Reassigment
[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

[tfr, tfrrs, tfrtic] = rsSTFTbase_RT(x, lowFreq, highFreq, alpha, hop, LL, h(1,:)', Dh(1,:)', 1);


%% Multitapering
ConceFT = [] ;

if MT > 1

	% Conceft
    ConceFT = zeros(size(tfrrs)) ;

    for ii = 1: MT
		rv = randn(1, dim) ; 
        rv = rv ./ norm(rv) ;
		rh = rv * h ; 
		rDh = rv * Dh ;

		[~, tfrrs, ~] = rsSTFTbase_RT(x, lowFreq, highFreq, alpha, hop, LL, rh', rDh', 1);

	 	ConceFT = ConceFT + tfrrs ;
    end

    ConceFT = ConceFT ./ MT ;

end

end
