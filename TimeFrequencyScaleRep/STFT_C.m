function [tfr, tfrtic] = STFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp)
% STFT_C  STFT using Hermite window
% Usage: 
% 	[tfr, tfrtic] = STFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp)
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
% 
% Output:
%   tfr: STFT
%   tfrtic: frequencies for which STFT is evaluated

[h, ~, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

[tfr, tfrtic] = STFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)');