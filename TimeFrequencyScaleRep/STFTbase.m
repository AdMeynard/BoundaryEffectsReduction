function [tfr, tfrtic] = STFTbase(x, lowFreq, highFreq, alpha, tDS, h)

[xrow,xcol] = size(x) ;
t = 1:length(x) ;
tLen = length(t(1:tDS:length(x))) ;

	% for tfr
N = length(-0.5+alpha:alpha:0.5) ;

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
	%% run STFT
tfr = zeros(fLen, tLen) ;
tfrtic = linspace(lowFreq, highFreq, fLen)' ;

for tidx = 1:tLen

    ti = t((tidx-1)*tDS+1); 
    tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices = rem(N+tau,N)+1;
    norm_h=norm(h(Lh+1+tau));

	tf0 = zeros(N, 1) ;
    tf0(indices) = x(ti+tau).*conj( h(Lh+1+tau)) /norm_h;
    tf0 = fft(tf0) ;

	tfr(:, tidx) = tf0(Lidx:Hidx) ;

end