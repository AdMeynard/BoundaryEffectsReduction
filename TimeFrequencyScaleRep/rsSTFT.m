function [tfr, tfrrs, tfrtic] = rsSTFT(x, lowFreq, highFreq, alpha, tDS, WinLen, dim, supp)
% RSSTFTBASE Reassigment of the STFT using Hermite window
% Usage: 
% 	[tfr, tfrrs, tfrtic] = rsSTFT(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp)
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
%   tfrrs: Reassigment of STFT
%   ConceFT: ConceFT RS of STFT
%   tfrtic: frequencies for which TFRs are evaluated


[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)
[tfr, tfrrs, tfrtic] = rsSTFTbase(x, lowFreq, highFreq, alpha, tDS, h(1,:)', Dh(1,:)', 1);

