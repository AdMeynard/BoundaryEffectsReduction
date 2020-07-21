function [forecastErr,CompTime, varargout] = SucForecast(xtot,fs,methods,HOP,N,extM,extK,extSEC,varargin)
% SUCFORECAST Run succesives forecasting and TFR on short subsignals of a large signal (for performance evaluation purpose)
% Usage:	[forecastErr,CompTime, OTD] = SucForecast(xtot,fs,methods,HOP,N,extM,extK,extSEC,tfr)
%
% Input:
%   xtot: large signal
%   fs: sampling frequency
%   methods (cell): list of forecasting methods to be implemented. See forecasting.m for the available methods   
%   HOP: subsampling rate for forecasting
%   N: segments length
%   extK: size of the dataset for forecasting
%   extM: lengths of the segments used for forecasting
%   extSEC: extension length in seconds
%   method: forcasting method. To be chosen between
%       'SigExt': Least square estimation
%       'edmd': Empirical dynamical Mode decompostion
%       'gpr': Gaussian process regression
%       'symmetrization': Symmetric extension
%   tfr (optionnal): time-frequency representations to be implemented. Available:
%       'conceFT': conceFT synchrosquezzing from STFT. Simultaneously provide the classic SST, and the STFT
%       'sst': synchrosquezzing from STFT. Simultaneously provide the STFT
%       'RS': reassigment from STFT.
% 
% Output:
%   forecastErr: forecasting mean square error for each segment
%   CompTime: computational time of the signal extension step, for each segment
%   OTD: optimal transport distance to the optimal TFR, for each segment

Ntot = length(xtot); % signal length
L = round(extSEC*fs) ;

if isempty(varargin)
    TFR = 'none' ;
else
    TFR = varargin{1} ;
    basicTF = varargin{2} ;
end

%% Forecastings

nleft = L+1 ;
nend = nleft+N+L-1 ;
ind = 1 ;

% figure;
while nend <= Ntot
    %% subsignal
    x = xtot(nleft:(nleft+N-1)) ;
    mu = mean(x);
    sigma = std(x);
    x = (x - mu) / sigma ;
    
    %% Extensions
    xxTRUE = ( xtot((nleft-L):nend) - mu ) / sigma ; % ground-truth extension
    
    if any( strcmp(methods, {'SigExt','lse'}) )
        method.name = 'SigExt' ;
        tic;
        xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; % Forecasted signal via LSE
        LSEtime(ind) = toc;
        forecastErr.LSE(ind) =  (1/(2*L)) * sqrt( sum( abs(xxLSE - xxTRUE).^2 ) ) ;
    end
    
    if any(strcmp(methods,'symmetrization'))
        method.name = 'symmetrization' ;
        tic;
        xxSYM = SigExtension(x,fs,HOP,extK,extM,extSEC,method) ; % Forecasted signal via LSE
        SYMtime(ind) = toc;
        forecastErr.SYM(ind) =  (1/(2*L)) * sqrt( sum( abs(xxSYM - xxTRUE).^2 ) ) ;
    end
    
    if any(strcmp(methods,'edmd'))
        method.name = 'edmd' ;
        method.param = 100 ;
        tic; 
        xxEDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
        EDMDtime(ind) = toc;
        forecastErr.EDMD(ind) = (1/(2*L)) * sqrt(sum(abs(xxEDMD - xxTRUE).^2)) ;
    end
    
    if any(strcmp(methods,'gpr'))
        method.name = 'gpr' ;
        tic; 
        xxGPR = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
        GPRtime(ind) = toc;
        forecastErr.GPR(ind) = (1/(2*L)) * sqrt(sum(abs(xxGPR - xxTRUE).^2)) ;
    end
    
    xxZP = [zeros(L,1); x; zeros(L,1)];
    forecastErr.ZP(ind) = (1/(2*L)) * sqrt(sum( abs(xxZP - xxTRUE).^2 )) ;

    %% Time-Frequency Representations
    if any(strcmp(TFR,'conceFT'))
        % SST parameters
        t = linspace(0, (N-1)/fs, N) ;
        tt = linspace(-L/fs, (N-1+L)/fs, N+2*L) ;
        tsst = t(1:basicTF.hop:end);
        tsstEXT = tt(1:basicTF.hop:end);
        
        n0 = find(tsstEXT>=min(tsst),1);
        nf = n0 + length(tsst)-1;

        % On the original long signal
        [VxxTRUE, sstTRUE, conceftTRUE, ~] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
        VxxTRUE = VxxTRUE(: , n0:nf )  ;
        sstTRUE = sstTRUE(: , n0:nf ) ;
        conceftTRUE = conceftTRUE(: , n0:nf ) ;
        % On the original short signal
        [Vx, sstS, conceftS, ~] = ConceFT_sqSTFT_C(double(x-mean(x)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
        OTD.stft.S(ind) = slicedOT(Vx, VxxTRUE) ;
        OTD.sst.S(ind) = slicedOT(sstS, sstTRUE) ;
        OTD.conceft.S(ind) = slicedOT(conceftS, conceftTRUE) ;
        % On the LS Vector estimated extended signal
        if any(strcmp(methods,{'SigExt','lse'}))
            [VxxLSE, sstLSE, conceftLSE, ~] = ConceFT_sqSTFT_C(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
            VxxLSE = VxxLSE(: , n0:nf );
            sstLSE = sstLSE(: , n0:nf );
            conceftLSE = conceftLSE(: , n0:nf );
            OTD.stft.S(ind) = slicedOT(VxxLSE, VxxTRUE) ;
            OTD.sst.LSE(ind) = slicedOT(sstLSE, sstTRUE) ;
            OTD.conceft.LSE(ind) = slicedOT(conceftLSE, conceftTRUE) ;
        end
        % On the SYM extended signal
        if any(strcmp(methods,'symmetrization'))
            [VxxSYM, sstSYM, conceftSYM, ~] = ConceFT_sqSTFT_C(double(xxSYM-mean(xxSYM)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
            VxxSYM = VxxSYM(:, n0:nf ) ;
            sstSYM = sstSYM(: , n0:nf ) ;
            conceftSYM = conceftSYM(: , n0:nf ) ;
            OTD.stft.SYM(ind) = slicedOT(VxxSYM, VxxTRUE) ;
            OTD.sst.SYM(ind) = slicedOT(sstSYM, sstTRUE) ;
            OTD.conceft.SYM(ind) = slicedOT(conceftSYM, conceftTRUE) ;
        end
        % On the EDMD Vector estimated extended signal
        if any(strcmp(methods,'edmd'))
            [VxxEDMD, sstEDMD, conceftEDMD, ~] = ConceFT_sqSTFT_C(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
            VxxEDMD = VxxEDMD(: , n0:nf ) ;
            sstEDMD = sstEDMD(: , n0:nf ) ;
            conceftEDMD = conceftEDMD(: , n0:nf ) ;
            OTD.stft.EDMD(ind) = slicedOT(VxxEDMD, VxxTRUE) ;
            OTD.sst.EDMD(ind) = slicedOT(sstEDMD, sstTRUE) ;
            OTD.conceft.EDMD(ind) = slicedOT(conceftEDMD, conceftTRUE) ;
        end
        % On the GPR Vector estimated extended signal
        if any(strcmp(methods,'gpr'))
            [VxxGPR, sstGPR, conceftGPR, ~] = ConceFT_sqSTFT_C(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
            VxxGPR = VxxGPR(: , n0:nf );
            sstGPR = sstGPR(: , n0:nf );
            OTD.stft.GPR(ind) = slicedOT(VxxGPR, VxxTRUE) ;
            conceftGPR = conceftGPR(: , n0:nf );
            OTD.sst.GPR(ind) = slicedOT(sstGPR, sstTRUE) ;
            OTD.conceft.GPR(ind) = slicedOT(conceftGPR, conceftTRUE) ;
        end
    end
    
    if any(strcmp(TFR,'sst'))
        % SST parameters
        t = linspace(0, (N-1)/fs, N) ;
        tt = linspace(-L/fs, (N-1+L)/fs, N+2*L) ;
        tsst = t(1:basicTF.hop:end);
        tsstEXT = tt(1:basicTF.hop:end);

        n0 = find(tsstEXT>=min(tsst),1);
        nf = n0 + length(tsst)-1;

        % On the original long signal
        [VxxTRUE, sstTRUE, ~] = sqSTFT_C(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
        sstTRUE = sstTRUE(: , n0:nf ) ; 
        VxxTRUE = VxxTRUE(: , n0:nf ) ;
        % On the original short signal
        [Vx, sstS, ~] = sqSTFT_C(double(x-mean(x)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
        OTD.sst.S(ind) = slicedOT(sstS, sstTRUE) ;
        OTD.stft.S(ind) = slicedOT(Vx, VxxTRUE) ;
        % On the LS Vector estimated extended signal
        if any(strcmp(methods,{'SigExt','lse'}))
            [VxxLSE, sstLSE, ~] = sqSTFT_C(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
            sstLSE = sstLSE(: , n0:nf );
            VxxLSE = VxxLSE(: , n0:nf );
            OTD.sst.LSE(ind) = slicedOT(sstLSE, sstTRUE) ;
            OTD.stft.LSE(ind) = slicedOT(VxxLSE, VxxTRUE) ;
        end
        % On the SYM Vector estimated extended signal
        if any(strcmp(methods,'symmetrization'))
            [VxxSYM, sstSYM, ~] = sqSTFT_C(double(xxSYM-mean(xxSYM)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
            sstSYM = sstSYM(: , n0:nf );
            VxxSYM = VxxSYM(: , n0:nf );
            OTD.sst.SYM(ind) = slicedOT(sstSYM, sstTRUE) ;
            OTD.stft.SYM(ind) = slicedOT(VxxSYM, VxxTRUE) ;
        end
        % On the EDMD Vector estimated extended signal
        if any(strcmp(methods,'edmd'))
            [VxxEDMD, sstEDMD, ~] = sqSTFT_C(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
            sstEDMD = sstEDMD(: , n0:nf );
            VxxEDMD = VxxEDMD(: , n0:nf );
            OTD.sst.EDMD(ind) = slicedOT(sstEDMD, sstTRUE) ;
            OTD.stft.EDMD(ind) = slicedOT(VxxEDMD, VxxTRUE) ;
        end
        % On the GPR Vector estimated extended signal
        if any(strcmp(methods,'gpr'))
            [VxxGPR, sstGPR, ~] = sqSTFT_C(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
            sstGPR = sstGPR(: , n0:nf );
            VxxGPR = VxxGPR(: , n0:nf );           
            OTD.sst.GPR(ind) = slicedOT(sstGPR, sstTRUE) ;
            OTD.stft.GPR(ind) = slicedOT(VxxGPR, VxxTRUE) ;
        end
    end
    
    if any(strcmp(TFR,'RS'))
        % SST parameters
        t = linspace(0, (N-1)/fs, N) ;
        tt = linspace(-L/fs, (N-1+L)/fs, N+2*L) ;
        trs = t(1:basicTF.hop:end) ;
        trsEXT = tt(1:basicTF.hop:end) ;
        n0 = find(trsEXT>=min(trs),1);
        nf = n0 + length(trs)-1;

        % On the original long signal
        [~, rsTRUE, ~] = rsSTFT(double(xxTRUE-mean(xxTRUE)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
        rsTRUE = rsTRUE(: , n0:nf ) ; 
        % On the original short signal
        [~, rsS, ~] = rsSTFT(double(x-mean(x)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
        OTD.rs.S(ind) = slicedOT(rsS, rsTRUE) ;
        % On the LS Vector estimated extended signal
        if any(strcmp(methods,{'SigExt','lse'}))
            [~, rsLSE, ~] = rsSTFT(double(xxLSE-mean(xxLSE)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
            rsLSE = rsLSE(: , n0:nf ) ;
            OTD.rs.LSE(ind) = slicedOT(rsLSE, rsTRUE) ;
        end
        % On the SYM Vector estimated extended signal
        if any(strcmp(methods,'symmetrization'))
            [~, rsSYM, ~] = rsSTFT(double(xxSYM-mean(xxSYM)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
            rsSYM = rsSYM(: , n0:nf ) ;
            OTD.rs.SYM(ind) = slicedOT(rsSYM, rsTRUE) ;
        end
        % On the EDMD Vector estimated extended signal
        if any(strcmp(methods,'edmd'))
            [~, rsEDMD, ~] = rsSTFT(double(xxEDMD-mean(xxEDMD)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
            rsEDMD = rsEDMD(: , n0:nf ) ;
            OTD.rs.EDMD(ind) = slicedOT(rsEDMD, rsTRUE) ;
        end
        % On the GPR Vector estimated extended signal
        if any(strcmp(methods,'gpr'))
            [~, rsGPR, ~] = rsSTFT(double(xxGPR-mean(xxGPR)), basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
            rsGPR = rsGPR(: , n0:nf ) ;
            OTD.rs.GPR(ind) = slicedOT(rsGPR, rsTRUE) ;
        end
    end
    
    nleft = nleft + N ;
    nend = nleft + N + L -1 ;
    ind = ind + 1 ;
end

if any(strcmp(methods,{'SigExt','lse'}))
    CompTime.LSE = mean(LSEtime) ;
end
if any(strcmp(methods,'symmetrization'))
    CompTime.SYM = mean(SYMtime) ;
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
