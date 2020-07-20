function [tfr, tfrsq, tfrtic] = sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, Second, Smooth)
% SQSTFT_C Synchrosqueezing of the STFT using Hermite window
% Usage: 
% 	[tfr, tfrsq, tfrtic] = sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, Second, Smooth)
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
%   Second: if 1, compute the second-order synchrosqueezing
%   Smooth: if 1, compute the smoothed version of the synchrosqueezing
% 
% Output:
%   tfr: STFT
%   tfrrs: Reassigment of STFT
%   ConceFT: ConceFT RS of STFT
%   tfrtic: frequencies for which TFRs are evaluated


[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

if ~Second
    [tfr, tfrsq, tfrtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', Smooth);
else
    [tfr, ~, tfrsq, tfrtic] = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', dwindow(Dh(1,:)'));
end
