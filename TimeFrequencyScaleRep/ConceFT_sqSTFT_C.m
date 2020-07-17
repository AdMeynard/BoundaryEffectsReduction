function [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Second, Smooth) ;

%
method = 'MEAN' ;
%method = 'MEDIAN' ;
N = length(x) ;

%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h, Dh, ~] = hermf(WinLen, dim, supp) ; % generate the window for short time Fourier transform (STFT)

if ~Second
    [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', Smooth);
else
    [tfr, tfrtic, ~, tfrsq, tfrsqtic] = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, h(1,:)', Dh(1,:)', dwindow(Dh(1,:)'));
end

%=======================================


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
			[~, ~, tfrsqX, tfrsqtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, rh', rDh', Smooth);
		else
			[~, ~, ~, tfrsqX, tfrsqtic] = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, rh', rDh', dwindow(rDh'));
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

	fprintf('\n') ;

end
