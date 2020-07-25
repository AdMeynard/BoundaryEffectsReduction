function [tfr, tfrtic] = STFT_RT(x, lowFreq, highFreq, alpha, hop, WinLen, LL, dim, supp)
% STFT_RT  STFT using Hermite window
% Usage: 
% 	[tfr, tfrtic] = STFT_RT(x, lowFreq, highFreq, alpha, hop, LL, WinLen, dim, supp)
%
% Input:
%   x: signal to be analized
%   lowFreq: relative lower frequency (0<lowFreq<Highfreq)
%   highFreq: relative higher frequency (lowFreq<Highfreq<0.5)
%   alpha: resolution in the frequency axis
%   hop: hop size (in samples)
%   LL: extension length
%   WinLen: Window length (in samples)
%   dim: order of the Hermite function used as window
%   supp: half time support of the Hermite function used as window
% 
% Output:
%   tfr: STFT
%   tfrtic: frequencies for which STFT is evaluated

[h, ~, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

[tfr, tfrtic] = STFTbase_RT(x, lowFreq, highFreq, alpha, hop, LL, h(1,:)');