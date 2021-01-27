function [TFRtot, ForecastTime, TFRtime] = BoundEffRed_RT(xtot,fs,forecastMethod,basicTF,VideoWriting,varargin)
% BOUNDEFFRED_RT Real-time implementation of BoundEffRed
% Usage:	[TFRtot, dt] = BoundEffRed_RT(xtot,fs,forecastMethod,basicTF,VideoWriting)
%
% Input:
%   xtot: signal to analyze
%   fs: sampling frequency
%   forecastMethod: forecasting method
%   basicTF: time-frequency representation and associated parameters
%   VideoWriting: When set to 1, the real-time TF represention image is recorded 
%
% Output:
%   TFRtot: Whole time-frequency representation
%   ForecastTime: forecastion iteration time
%   TFRtime: TF representation update time

if VideoWriting==1
    if isempty(varargin)
        error('The video name must be specified')
    else
        VideoName = ['../../Results/' varargin{1}] ;
    end
end


%% Parameters
Ntot = length(xtot) ;

HOP = 1 ;
L = round((basicTF.win-1)/2) ; % extension length
extM = round(1.5*L) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

Nt = extK + extM + 1 ; % segments length
LL = ceil(L/basicTF.hop); % extension in TF domain

freqs = fs * (basicTF.fmin:basicTF.df:basicTF.fmax) ;

%% Initialization
n0 = 1 ;
n1 = n0+Nt-1 ;

t = (n0:n1)/fs ;
x = xtot(n0:n1) ;

switch basicTF.representation
    case 'SST'
        [~, TFRcurrent, ~] = sqSTFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    case 'conceFT'
        [~, ~, TFRcurrent, ~] = ConceFT_sqSTFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
    case 'STFT'
        [TFRcurrent, ~] = STFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    case 'RS'
        [~, TFRcurrent, ~] = rsSTFT(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
end

TFRtot = TFRcurrent(:,end-LL) ;
        
figure;
imagesc(t,freqs,log1p(abs(TFRcurrent)/5e0));% colormap(1-gray);
axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
drawnow;

if VideoWriting
    myVideo = VideoWriter(VideoName);
    myVideo.FrameRate = 10 ;
    open(myVideo)
end

k = 1 ;
while n1<Ntot
    t = (n0:n1)/fs ;
    x = xtot(n0:n1) ;

    %% Forecasting
    tic;
    xext = forecasting(x,L,HOP,extK,extM,forecastMethod);
    xx = [x; xext]; % size Nt + L
    ForecastTime(k) = toc;
    
    %% Update TF Representation
    tic;
    xxPREC = xx((end-round(2.5*L)):end) ; % Part of the signal to be used to update TFR
    
    switch basicTF.representation
        case 'SST'
            [~, TFRxx, ~] = sqSTFT_RT(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 0, 0) ;
        case 'conceFT'
            [~, ~, TFRxx, ~] = ConceFT_sqSTFT_RT(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1, 0, 0) ;
        case 'STFT'
            [TFRxx, ~] = STFT_RT(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10) ;
        case 'RS'
            [~, TFRxx, ~, ~] = ConceFT_rsSTFT_RT(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, LL, 1, 10, 1) ;
    end

    TFRcurrent = [TFRcurrent(:,2:(end-LL)) TFRxx] ; % sliding sst
    
    TFRtime(k) = toc ;
    
    %% Display
    
    imagesc(t,freqs,log1p(abs(TFRcurrent)/5e0)); colormap(1-gray); 
    axis xy; % xlabel('Time (s)'); ylabel('Frequency (Hz)');
    drawnow;
    
    n0 = n0 + basicTF.hop ;
    n1 = n0 + Nt - 1 ;
    
    TFRtot = [TFRtot TFRxx(:,end)] ; % whole SST (non real-time)
    k = k+1 ;
    
    if VideoWriting
        frame = getframe(gcf) ;
        writeVideo(myVideo, frame);
    end
    
end

TFRtot = [TFRtot TFRxx(:,(end-LL+1):end)] ;

xlabel('Time (s)'); ylabel('Frequency (Hz)');

if VideoWriting
    close(myVideo) ;
end
