function [forecastErr,CompTime, varargout] = SucForecast(xtot,fs,methods,HOP,N,extM,extK,extSEC,varargin)
% Run succesives forecasting and SST on short subsignals of a large signal

Ntot = length(xtot); % signal length
L = round(extSEC*fs) ;

if isempty(varargin)
    TFR = 'none' ;
else
    TFR = varargin{1} ;
    basicTF = varargin{2} ;
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
    xxTRUE = ( xtot((k*N+1-L):nend) - mu ) / sigma ; % ground-truth extension
    
    if any(strcmp(methods,'lse'))
        method.name = 'lse' ;
        tic;
        xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; % Forecasted signal via LSE
        LSEtime(ind) = toc;
        forecastErr.LSE(ind) = sqrt( (1/(2*L)) * sum( abs(xxLSE - xxTRUE).^2 ) ) ;
    end
    
    if any(strcmp(methods,'edmd'))
        method.name = 'edmd' ;
        method.param = 100 ;
        tic; 
        xxEDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
        EDMDtime(ind) = toc;
        forecastErr.EDMD(ind) = (1/(2*extSEC*fs)) * sqrt(sum(abs(xxEDMD - xxTRUE).^2)) ;
    end
    
    if any(strcmp(methods,'gpr'))
        method.name = 'gpr' ;
        tic; 
        xxGPR = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
        GPRtime(ind) = toc;
        forecastErr.GPR(ind) = (1/(2*extSEC*fs)) * sqrt(sum(abs(xxGPR - xxTRUE).^2)) ;
    end
    
    xxZP = [zeros(L,1); x; zeros(L,1)];
    forecastErr.ZP(ind) = sqrt( (1/(2*L)) * sum( abs(xxZP - xxTRUE).^2 ) ) ;

    %% Time-Frequency Representations
    if any(strcmp(TFR,'conceFT'))
        % SST parameters
        t = linspace(0, (N-1)/fs, N) ;
        tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
        tsst = t(1:basicTF.hop:end);
        tsstEXT = tt(1:basicTF.hop:end);

        % On the original long signal
        [~, ~, sstTRUE, conceftTRUE, ~] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        sstTRUE = sstTRUE(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
        conceftTRUE = conceftTRUE(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
        % On the original short signal
        [~, ~, sstS, conceftS, ~] = ConceFT_sqSTFT_C(double(x-mean(x)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        OTD.sst.S(ind) = slicedOT(sstS, sstTRUE) ;
        OTD.conceft.S(ind) = slicedOT(conceftS, conceftTRUE) ;
        % On the LS Vector estimated extended signal
        if any(strcmp(methods,'lse'))
            [~, ~, sstLSE, conceftLSE, ~] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstLSE = sstLSE(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            conceftLSE = conceftLSE(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            OTD.sst.LSE(ind) = slicedOT(sstLSE, sstTRUE) ;
            OTD.conceft.LSE(ind) = slicedOT(conceftLSE, conceftTRUE) ;
        end
        % On the EDMD Vector estimated extended signal
        if any(strcmp(methods,'edmd'))
            [~, ~, sstEDMD, conceftEDMD, ~] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstEDMD = sstEDMD(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            conceftEDMD = conceftEDMD(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            OTD.sst.EDMD(ind) = slicedOT(sstEDMD, sstTRUE) ;
            OTD.conceft.EDMD(ind) = slicedOT(conceftEDMD, conceftTRUE) ;
        end
        % On the GPR Vector estimated extended signal
        if any(strcmp(methods,'gpr'))
            [~, ~, sstGPR, ~, ~] = ConceFT_sqSTFT_C(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstGPR = sstGPR(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            OTD.sst.GPR(ind) = slicedOT(sstGPR, sstTRUE) ;
            conceftGPR = conceftGPR(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
            OTD.conceft.GPR(ind) = slicedOT(conceftGPR, conceftTRUE) ;
        end
    end
    
    if any(strcmp(TFR,'sstSTFT'))
        % SST parameters
        t = linspace(0, (N-1)/fs, N) ;
        tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
        tsst = t(1:basicTF.hop:end);
        tsstEXT = tt(1:basicTF.hop:end);

        n0 = find(tsstEXT>=min(tsst),1);
        nf = n0 + length(tsst)-1;

        % On the original long signal
        [VxxTRUE, ~, sstTRUE, ~, ~] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        sstTRUE = sstTRUE(: , n0:nf ); 
        VxxTRUE = VxxTRUE(: , n0:nf );
        % On the original short signal
        [Vx, ~, sstS, ~] = ConceFT_sqSTFT_C(double(x-mean(x)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        OTD.sst.S(ind) = slicedOT(sstS, sstTRUE) ;
        OTD.stft.S(ind) = slicedOT(Vx, VxxTRUE) ;
        % On the LS Vector estimated extended signal
        if any(strcmp(methods,'lse'))
            [VxxLSE, ~, sstLSE, ~, ~] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstLSE = sstLSE(: , n0:nf );
            VxxLSE = VxxLSE(: , n0:nf );
            OTD.sst.LSE(ind) = slicedOT(sstLSE, sstTRUE) ;
            OTD.stft.LSE(ind) = slicedOT(VxxLSE, VxxTRUE) ;
        end
        % On the EDMD Vector estimated extended signal
        if any(strcmp(methods,'edmd'))
            [VxxEDMD, ~, sstEDMD, ~, ~] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstEDMD = sstEDMD(: , n0:nf );
            VxxEDMD = VxxEDMD(: , n0:nf );
            OTD.sst.EDMD(ind) = slicedOT(sstEDMD, sstTRUE) ;
            OTD.stft.EDMD(ind) = slicedOT(VxxEDMD, VxxTRUE) ;
        end
        % On the GPR Vector estimated extended signal
        if any(strcmp(methods,'gpr'))
            [VxxGPR, ~, sstGPR, ~, ~] = ConceFT_sqSTFT_C(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
            sstGPR = sstGPR(: , n0:nf );
            VxxGPR = VxxGPR(: , n0:nf );           
            OTD.sst.GPR(ind) = slicedOT(sstGPR, sstTRUE) ;
            OTD.stft.GPR(ind) = slicedOT(VxxGPR, VxxTRUE) ;
        end
    end
    
    if any(strcmp(TFR,'RS'))
        % SST parameters
        t = linspace(0, (N-1)/fs, N) ;
        tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
        trs = t;
        trsEXT = tt ;
        n0 = find(trsEXT>=min(trs),1);
        nf = n0 + length(trs)-1;

        % On the original long signal
        [~, ~, rsTRUE, ~, ~] = ConceFT_rsSTFT(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
        rsTRUE = rsTRUE(: , n0:nf ) ; 
        % On the original short signal
        [~, ~, rsS, ~] = ConceFT_rsSTFT(double(x-mean(x)), basicTF.fmin, basicTF.fmax, 1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
        OTD.rs.S(ind) = slicedOT(rsS, rsTRUE) ;
        % On the LS Vector estimated extended signal
        if any(strcmp(methods,'lse'))
            [~, ~, rsLSE, ~, ~] = ConceFT_rsSTFT(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsLSE = rsLSE(: , n0:nf ) ;
            OTD.rs.LSE(ind) = slicedOT(rsLSE, rsTRUE) ;
        end
        % On the EDMD Vector estimated extended signal
        if any(strcmp(methods,'edmd'))
            [~, ~, rsEDMD, ~, ~] = ConceFT_rsSTFT(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsEDMD = rsEDMD(: , n0:nf ) ;
            OTD.rs.EDMD(ind) = slicedOT(rsEDMD, rsTRUE) ;
        end
        % On the GPR Vector estimated extended signal
        if any(strcmp(methods,'gpr'))
            [~, ~, rsGPR, ~, ~] = ConceFT_rsSTFT(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax,1e-5, basicTF.hop, basicTF.win, 1, 10, 1) ;
            rsGPR = rsGPR(: , n0:nf ) ;
            OTD.rs.GPR(ind) = slicedOT(rsGPR, rsTRUE) ;
        end
    end
    
    k = k+1 ;
    nend = (k+1)*N + extSEC*fs ;
    ind = ind + 1 ;
end

if any(strcmp(methods,'lse'))
    CompTime.LSE = mean(LSEtime) ;
end
if any(strcmp(methods,'edmd'))
    CompTime.EDMD = mean(EDMDtime) ;
end
if any(strcmp(methods,'gpr'))
    CompTime.GPR = mean(GPRtime) ;
end

if ~strcmpi(TFR,'none') 
    varargout{1} = OTD ;
else
    varargout{1} = {};
end
