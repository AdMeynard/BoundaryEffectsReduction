function [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, Second, Smooth)

[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

if ~Second
    [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', Smooth);
else
    [tfr, tfrtic, ~, tfrsq, tfrsqtic] = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', dwindow(Dh(1,:)'));
end
