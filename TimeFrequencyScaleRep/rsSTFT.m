function [tfr, tfrtic, tfrrs, tfrrstic] = rsSTFT(x, lowFreq, highFreq, alpha, tDS, WinLen, dim, supp)
%
% Usage: 
% 	[tfr, tfrtic, tfrrs, tfrrstic] = rSTFT(x, lowFreq, highFreq, alpha, WinLen, dim, supp)
%
% alpha: resolution in the frequency axis
% WinLen, dim, supp: for hermf.m
%
% Example:
% 	[tfr, tfrtic, tfrrs, tfrrstic] = rsSTFT(y, 0,0.5, 0.0002, 121, 4, 6, 10);


[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)
[tfr, tfrtic, tfrrs, tfrrstic] = rsSTFTbase(x, lowFreq, highFreq, alpha, tDS, h(1,:)', Dh(1,:)', 1);

