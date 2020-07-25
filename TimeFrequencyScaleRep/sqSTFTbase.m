function [tfr, tfrsq, tfrtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, tDS, h, Dh, Smooth)
% SQSTFTBASE Synchrosqueezing of the STFT
% Usage: 
% 	[tfr, tfrrs, tfrtic] = sqSTFTbase(x, lowFreq, highFreq, alpha, tDS, h, Dh, Smooth)
%
% Input:
%	x     : analysed signal.
%	[lowFreq, highFreq] : frequency range \subset [0 0.5]. For the sake of computational efficiency.
%	alpha : the resolution in the frequency axis
%	tDS   : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
%	h     : frequency smoothing window, H(0) being forced to 1
%   Dh    : differentiation of H
%   Smooth: if 1, compute a smoothed version of the synchrosqueezing
% 
% Output:
%	tfr   : STFT
%	tfrrs : synchrosqueezed STFT
%   tfrtic: frequencies for which TFRs are evaluated
%
%	F. Auger, May-July 1994, July 1995.  
%  modifed from tfrrsp.m, by Hau-tieng Wu, 2013 
%	Copyright (c) 1996 by CNRS (France).
%

[xrow,xcol] = size(x) ;
t = 1:length(x) ;
tLen = length(t(1:tDS:length(x))) ;

	% for tfr
N = length(-0.5+alpha:alpha:0.5) ;

	% for tfrsq
Lidx = round( (N/2)*(lowFreq/0.5) ) + 1 ; 
Hidx = round( (N/2)*(highFreq/0.5) ) ; 
fLen = Hidx - Lidx + 1 ;

%====================================================================
	%% check input signals
if (xcol~=1)
    error('X must have only one column');
elseif highFreq > 0.5
    error('TopFreq must be a value in [0, 0.5]');
elseif (tDS < 1) || (rem(tDS,1)) 
    error('tDS must be an integer value >= 1');
end

[hrow,hcol] = size(h); Lh = (hrow-1)/2; 
if (hcol~=1)||(rem(hrow,2)==0)
    error('H must be a smoothing window with odd length');
end


%====================================================================
	%% run STFT and reassignment rule
tfr = zeros(fLen, tLen) ;
tfrsq = zeros(fLen, tLen) ; 
tfrtic = linspace(lowFreq, highFreq, fLen)' ;

Ex = mean(abs(x).^2);
Threshold = 1.0e-8*Ex;  % originally it was 1e-6*Ex

Mid = round(length(tfrtic)/2) ;
Delta = 20*(tfrtic(2)-tfrtic(1)).^2 ;
weight = exp(-(tfrtic(Mid-10:Mid+10)-tfrtic(Mid)).^2/Delta) ;
weight = weight ./ sum(weight) ;
weightIDX = (Mid-10:Mid+10) - Mid ;

for tidx = 1:tLen

    ti = t((tidx-1)*tDS+1); 
    tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;
    norm_h=norm(h(Lh+1+tau));

	tf0 = zeros(N, 1) ; tf1 = zeros(N, 1) ;
    tf0(indices) = x(ti+tau).*conj( h(Lh+1+tau)) /norm_h;
    tf1(indices) = x(ti+tau).*conj(Dh(Lh+1+tau)) /norm_h;
    tf0 = fft(tf0) ;
    tf1 = fft(tf1) ;

		% get the first order omega
	omega = zeros(size(tf1)) ;
	avoid_warn = find(tf0~=0);
	omega(avoid_warn) = round( imag(N*tf1(avoid_warn)./tf0(avoid_warn)/(2.0*pi)) ) ;

	sst = zeros(fLen,1) ;

    for jcol = 1: N/2
        if abs(tf0(jcol)) > Threshold

   	    	jcolhat = jcol - omega(jcol) ;

            if ((jcolhat <= Hidx) && (jcolhat >= Lidx))
                
                if Smooth
					IDXb = find((jcolhat-Lidx+1+weightIDX <= Hidx) & (jcolhat-Lidx+1+weightIDX >= Lidx)) ;
					IDXa = jcolhat-Lidx+1+weightIDX(IDXb) ; 
                   	sst(IDXa) = sst(IDXa) + tf0(jcol)*weight(IDXb) ;
                else
	   				sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) + tf0(jcol) ;
                end
                
            end
        end
    end

	tfr(:, tidx) = tf0(Lidx:Hidx) ;
	tfrsq(:, tidx) = sst ;

end
