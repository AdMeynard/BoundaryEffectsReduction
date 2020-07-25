function [tfr, tfrsq, ConceFT, tfrtic] = ConceFT_sqSTFT_RT(x, lowFreq, highFreq, alpha, hop, WinLen, LL, dim, supp, MT, Second, Smooth)
% CONCEFT_SQSTFT_RT ConceFT Synchrosqueezing of the STFT
% Usage: 
% 	[tfr, tfrsq, ConceFT, tfrtic] = ConceFT_sqSTFT_RT(x, lowFreq, highFreq, alpha, hop, WinLen, LL, dim, supp, MT, Second, Smooth)
%
% Input:
%   x: signal to be analized
%   lowFreq: relative lower frequency (0<lowFreq<Highfreq)
%   highFreq: relative higher frequency (lowFreq<Highfreq<0.5)
%   alpha: resolution in the frequency axis
%   hop: hop size (in samples)
%   WinLen: Window length (in samples)
%   LL: extension length (in samples)
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

%% Ordinary SST
[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

if ~Second
    [tfr, tfrsq, tfrtic] = sqSTFTbase_RT(x, lowFreq, highFreq, alpha, hop, LL, h(1,:)', Dh(1,:)', Smooth);
else
    [tfr, ~, tfrsq, tfrtic] = sqSTFTbase2nd_RT(x, lowFreq, highFreq, alpha, hop, LL, h(1,:)', Dh(1,:)', dwindow(Dh(1,:)'));
end

%% Multitapering

ConceFT = abs(tfrsq) ;

if MT > 1

	% Conceft
    for ii = 1: MT
		rv = randn(1, dim) + 1i*randn(1, dim) ; 
        rv = rv ./ norm(rv) ;
		rh = rv * h ; 
		rDh = rv * Dh ;

		if ~Second
			[~, tfrsqX, ~] = sqSTFTbase_RT(x, lowFreq, highFreq, alpha, hop, LL, rh', rDh', Smooth);
		else
			[~, ~, tfrsqX, ~] = sqSTFTbase2nd_RT(x, lowFreq, highFreq, alpha, hop, LL, rh', rDh', dwindow(rDh'));
		end

	 	ConceFT = ConceFT + abs(tfrsqX) ;
    end

    ConceFT = ConceFT ./ (MT+1) ;

end