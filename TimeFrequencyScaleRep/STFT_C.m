function [tfr, tfrtic] = STFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp)

% generate the window for short time Fourier transform (STFT)
[h, ~, ~] = hermf(WinLen, dim, supp) ;
[tfr, tfrtic] = STFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)');

end
