function [forecastErr, varargout] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,varargin)
% Run succesives forecasting and SST on short subsignals of a large signal

Ntot = length(xtot); % signal length
L = round(extSEC*fs) ;

if isempty(varargin)
    basicTF.meth = 'none' ;
else
    basicTF = varargin{1} ;
end

%% Forecastings

k = ceil( (extSEC*fs)/(N+1) ) ;
nend = (k+1)*N + L ;
ind = 1 ;

% figure;
while nend <= Ntot
    % subsignal
    x = xtot((k*N+1):((k+1)*N)) ;
    mu = mean(x);
    sigma = std(x);
    x = (x - mu) / sigma ;
    
    % Extensions
    method.name = 'lse' ;
    xxLSEV = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; % Forecasted signal via LSE
%     xxDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,'dmd') ; % Forecasted signal via LSE
    xxZP = [zeros(L,1); x; zeros(L,1)];

    xxTRUE = ( xtot((k*N+1-L):nend) - mu ) / sigma ;
    
    
%     forecastErr.LSE(ind) =  sqrt( (1/(2*L)) * sum( abs(xxLSE - xxTRUE).^2 ) ) ;
    forecastErr.LSEV(ind) = sqrt( (1/(2*L)) * sum( abs(xxLSEV - xxTRUE).^2 ) ) ;
    forecastErr.ZP(ind) = sqrt( (1/(2*L)) * sum( abs(xxZP - xxTRUE).^2 ) ) ;
%     forecastErr.DMD(ind) = (1/(2*extSEC*fs)) * sqrt(sum(abs(xxDMD - xxTRUE).^2)) ;


    %% SST
    switch basicTF.meth
        case 'sst'
            % SST parameters
            t = linspace(0, (N-1)/fs, N) ;
            tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
            tsst = t(1:basicTF.hop:end);
            tsstEXT = tt(1:basicTF.hop:end);

            % On the original short signal
            [~, ~, sstS, ~] = ConceFT_sqSTFT_C(double(x-mean(x)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            % On the original long signal
            [~, ~, sstTRUE, ~, ~] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstTRUE = sstTRUE(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) ); 
             % On the LS Vector estimated extended signal 
            [~, ~, sstLSEV, ~, ~] = ConceFT_sqSTFT_C(double(xxLSEV-mean(xxLSEV)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstLSEV = sstLSEV(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            
            OTDS(ind) = slicedOT(sstS, sstTRUE) ;
            OTDLSE(ind) = slicedOT(sstLSE, sstTRUE) ;
            
        case 'sstSTFT'
            % SST parameters
            t = linspace(0, (N-1)/fs, N) ;
            tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
            tsst = t(1:basicTF.hop:end);
            tsstEXT = tt(1:basicTF.hop:end);
            
            n0 = find(tsstEXT>=min(tsst),1);
            nf = n0 + length(tsst)-1;

            % On the original short signal
            [Vx, ~, sstS, ~] = ConceFT_sqSTFT_C(double(x-mean(x)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            % On the original long signal
            [VxxTRUE, ~, sstTRUE, ~, ~] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstTRUE = sstTRUE(: , n0:nf ); 
            VxxTRUE = VxxTRUE(: , n0:nf );
             % On the LS Vector estimated extended signal 
            [VxxLSEV, ~, sstLSEV, ~, ~] = ConceFT_sqSTFT_C(double(xxLSEV-mean(xxLSEV)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstLSEV = sstLSEV(: , n0:nf );
            VxxLSEV = VxxLSEV(: , n0:nf );

            sstOTDS(ind) = slicedOT(sstS, sstTRUE) ;
            sstOTDLSE(ind) = slicedOT(sstLSEV, sstTRUE) ;
            
            stftOTDS(ind) = slicedOT(Vx, VxxTRUE) ;
            stftOTDLSE(ind) = slicedOT(VxxLSEV, VxxTRUE) ;
            
        case 'RS'
            % SST parameters
            t = linspace(0, (N-1)/fs, N) ;
            tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
            trs = t;
            trsEXT = tt ;
            n0 = find(trsEXT>=min(trs),1);
            nf = n0 + length(trs)-1;

            % On the original short signal
            [~, ~, rsS, ~] = ConceFT_rsSTFT(double(x-mean(x)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            % On the original long signal
            [~, ~, rsTRUE, ~, ~] = ConceFT_rsSTFT(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsTRUE = rsTRUE(: , n0:nf ); 
             % On the LS Vector estimated extended signal 
            [~, ~, rsLSE, ~, ~] = ConceFT_rsSTFT(double(xxLSEV-mean(xxLSEV)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsLSE = rsLSE(: , n0:nf );


            OTDS(ind) = slicedOT(rsS, rsTRUE) ;
            OTDLSE(ind) = slicedOT(rsLSE, rsTRUE) ;
        otherwise
            % nothing to do
    end
    
    k = k+1 ;
    nend = (k+1)*N + extSEC*fs ;
    ind = ind + 1 ;
end

switch basicTF.meth
    case {'sst','RS'}
        varargout{1} = OTDS ;
        varargout{2} = OTDLSE ;
    case 'sstSTFT'
        varargout{1} = sstOTDS ;
        varargout{2} = sstOTDLSE ;
        varargout{3} = stftOTDS ;
        varargout{4} = stftOTDLSE ;  
    case 'none'
        varargout{1} = {};
end
