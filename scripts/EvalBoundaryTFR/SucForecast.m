function [forecastErr,CompTime, varargout] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,varargin)
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
    %% subsignal
    x = xtot((k*N+1):((k+1)*N)) ;
    mu = mean(x);
    sigma = std(x);
    x = (x - mu) / sigma ;
    
    %% Extensions
    method.name = 'lse' ;
    tic;
    xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; % Forecasted signal via LSE
    LSEtime(ind) = toc;
    
    method.name = 'edmd' ;
    method.param = 100 ;
    tic; 
    xxEDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
    EDMDtime(ind) = toc;
    
    method.name = 'gpr' ;
    tic; 
    xxGPR = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
    GPRtime(ind) = toc;
    
    xxZP = [zeros(L,1); x; zeros(L,1)];

    xxTRUE = ( xtot((k*N+1-L):nend) - mu ) / sigma ;
    
    
%     forecastErr.LSE(ind) =  sqrt( (1/(2*L)) * sum( abs(xxLSE - xxTRUE).^2 ) ) ;
    forecastErr.LSE(ind) = sqrt( (1/(2*L)) * sum( abs(xxLSE - xxTRUE).^2 ) ) ;
    forecastErr.ZP(ind) = sqrt( (1/(2*L)) * sum( abs(xxZP - xxTRUE).^2 ) ) ;
    forecastErr.EDMD(ind) = (1/(2*extSEC*fs)) * sqrt(sum(abs(xxEDMD - xxTRUE).^2)) ;
    forecastErr.GPR(ind) = (1/(2*extSEC*fs)) * sqrt(sum(abs(xxGPR - xxTRUE).^2)) ;


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
            [~, ~, sstLSE, ~, ~] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstLSE = sstLSE(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            % On the LS Vector estimated extended signal 
            [~, ~, sstEDMD, ~, ~] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstEDMD = sstEDMD(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            % On the LS Vector estimated extended signal 
            [~, ~, sstGPR, ~, ~] = ConceFT_sqSTFT_C(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstGPR = sstGPR(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            
            OTD.sst.S(ind) = slicedOT(sstS, sstTRUE) ;
            OTD.sst.LSE(ind) = slicedOT(sstLSE, sstTRUE) ;
            OTD.sst.EDMD(ind) = slicedOT(sstEDMD, sstTRUE) ;
            OTD.sst.GPR(ind) = slicedOT(sstGPR, sstTRUE) ;
            
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
            [VxxLSE, ~, sstLSE, ~, ~] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstLSE = sstLSE(: , n0:nf );
            VxxLSE = VxxLSE(: , n0:nf );
            % On the EDMD Vector estimated extended signal 
            [VxxEDMD, ~, sstEDMD, ~, ~] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstEDMD = sstEDMD(: , n0:nf );
            VxxEDMD = VxxEDMD(: , n0:nf );
            % On the LS Vector estimated extended signal 
            [VxxGPR, ~, sstGPR, ~, ~] = ConceFT_sqSTFT_C(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstGPR = sstGPR(: , n0:nf );
            VxxGPR = VxxGPR(: , n0:nf );

            OTD.sst.S(ind) = slicedOT(sstS, sstTRUE) ;
            OTD.sst.LSE(ind) = slicedOT(sstLSE, sstTRUE) ;
            OTD.sst.EDMD(ind) = slicedOT(sstEDMD, sstTRUE) ;
            OTD.sst.GPR(ind) = slicedOT(sstGPR, sstTRUE) ;
            
            OTD.stft.S(ind) = slicedOT(Vx, VxxTRUE) ;
            OTD.stft.LSE(ind) = slicedOT(VxxLSE, VxxTRUE) ;
            OTD.stft.EDMD(ind) = slicedOT(VxxEDMD, VxxTRUE) ;
            OTD.stft.GPR(ind) = slicedOT(VxxGPR, VxxTRUE) ;
            
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
            [~, ~, rsLSE, ~, ~] = ConceFT_rsSTFT(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsLSE = rsLSE(: , n0:nf );
            % On the EDMD Vector estimated extended signal 
            [~, ~, rsEDMD, ~, ~] = ConceFT_rsSTFT(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsEDMD = rsEDMD(: , n0:nf );
            % On the GPR Vector estimated extended signal 
            [~, ~, rsGPR, ~, ~] = ConceFT_rsSTFT(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsGPR = rsGPR(: , n0:nf );


            OTD.rs.S(ind) = slicedOT(rsS, rsTRUE) ;
            OTD.rs.LSE(ind) = slicedOT(rsLSE, rsTRUE) ;
            OTD.rs.EDMD(ind) = slicedOT(rsEDMD, rsTRUE) ;
            OTD.rs.GPR(ind) = slicedOT(rsEDMD, rsTRUE) ;
        otherwise
            % nothing to do
    end
    
    k = k+1 ;
    nend = (k+1)*N + extSEC*fs ;
    ind = ind + 1 ;
end

CompTime.LSE = mean(LSEtime) ;
CompTime.EDMD = mean(EDMDtime) ;
CompTime.GPR = mean(GPRtime) ;

if ~strcmpi(basicTF.meth,'none') 
    varargout{1} = OTD ;
else
    varargout{1} = {};
end
