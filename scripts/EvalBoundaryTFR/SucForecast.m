function [forecastErr, varargout] = SucForecast(xtot,fs,HOP,N,extM,extK,extSEC,varargin)
% Run succesives forecasting and SST on short subsignals of a large signal

Ntot = length(xtot); % signal length
L = round(extSEC*fs) ;

if ~isempty(varargin)
    isSST = strcmpi( varargin{1}, 'SST') ; % check is SST computions are required
else
    isSST = 0 ;
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
%     xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,'lse') ; % Forecasted signal via LSE
    xxLSEV = SigExtension(x,fs,HOP,extK,extM,extSEC,'lseV') ; % Forecasted signal via LSE
%     xxDMD = SigExtension(x,fs,HOP,extK,extM,extSEC,'dmd') ; % Forecasted signal via LSE
    xxZP = [zeros(L,1); x; zeros(L,1)];

    xxTRUE = ( xtot((k*N+1-L):nend) - mu ) / sigma ;
    
    
%     forecastErr.LSE(ind) =  sqrt( (1/(2*L)) * sum( abs(xxLSE - xxTRUE).^2 ) ) ;
    forecastErr.LSEV(ind) = sqrt( (1/(2*L)) * sum( abs(xxLSEV - xxTRUE).^2 ) ) ;
    forecastErr.ZP(ind) = sqrt( (1/(2*L)) * sum( abs(xxZP - xxTRUE).^2 ) ) ;
%     forecastErr.DMD(ind) = (1/(2*extSEC*fs)) * sqrt(sum(abs(xxDMD - xxTRUE).^2)) ;

%     figure;
%     plot(tt,xxLSEV,tt,xxTRUE,'--',t,x,'linewidth',2); grid on;
%     legend('Estimated Extended signal','Ground truth Extended signal','Original signal'); 
%     xlabel('Time (s)'); ylabel('Signals'); title('Time series'); drawnow;

    %% SST
    if isSST
        % SST parameters
        basicTF.hop = 10;
        basicTF.win = 1501;
        t = linspace(0, (N-1)/fs, N) ;
        tt = linspace(-extSEC, (N-1)/fs+extSEC, N+2*extSEC*fs) ;
        tsst = t(1:basicTF.hop:end);
        tsstEXT = tt(1:basicTF.hop:end);
        
        % On the original short signal
        [~, ~, tfrsq3, ~, ~] = ConceFT_sqSTFT_C(double(x-mean(x)), 0, 0.01,...
                    1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        % On the original long signal
        [~, ~, tfrsq3TRUE, ~, ~] = ConceFT_sqSTFT_C(double(xxTRUE-mean(xxTRUE)), 0, 0.01,...
                    1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        tfrsq3TRUE = tfrsq3TRUE(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) ); 
         % On the LS Vector estimated extended signal 
        [~, ~, tfrsq3LSEV, ~, ~] = ConceFT_sqSTFT_C(double(xxLSEV-mean(xxLSEV)), 0, 0.01,...
                    1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        tfrsq3LSEV = tfrsq3LSEV(: , (tsstEXT>=min(tsst) & tsstEXT<=max(tsst)) );
        % On the zero-padded extended signal 
        [~, ~, tfrsq3ZP, ~, ~] = ConceFT_sqSTFT_C(double(xxZP-mean(xxZP)), 0, 0.01,...
                1e-5, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;


        sstS{ind} = tfrsq3 ;
        sstTRUE{ind} = tfrsq3TRUE ;
        sstLSE{ind} = tfrsq3LSEV ;
        sstZP{ind} = tfrsq3ZP ;
    end
    
    k = k+1 ;
    nend = (k+1)*N + extSEC*fs ;
    ind = ind + 1 ;
end

if isSST
    varargout{1} = sstS ;
    varargout{2} = sstTRUE ;
    varargout{3} = sstLSE ;
    varargout{4} = sstZP ;
else
    varargout{1} = {};
end