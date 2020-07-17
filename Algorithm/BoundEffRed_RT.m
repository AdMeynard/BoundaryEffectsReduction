function [TFRtot, dt] = BoundEffRed_RT(xtot,fs,forecastMethod,basicTF,VideoWriting)
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
%   dt: iteration time


%% Parameters
Ntot = length(xtot) ;

HOP = 1 ;
L = round((basicTF.win-1)/2) ; % extension length
extM = round(1.5*L) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

Nt = extK + extM + 1 ; % segments length
Lo = round(1.5*L/basicTF.hop); % overlap in TF domain
LL = round(L/basicTF.hop); % extension in TF domain

freqs = fs * (basicTF.fmin:basicTF.df:basicTF.fmax) ;

%% Initialization
n0 = 1 ;
n1 = n0+Nt-1 ;

t = (n0:n1)/fs ;
x = xtot(n0:n1) ;

switch basicTF.representation
    case 'SST'
        [~, ~, TFRcurrent, ~] = sqSTFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
    case 'conceFT'
        [~, ~, ~, TFRcurrent, ~] = ConceFT_sqSTFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, basicTF.MT, 0, 0) ;
    case 'STFT'
        [TFRcurrent, ~] = STFT_C(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
    case 'RS'
        [~, ~, TFRcurrent, ~, ~] = ConceFT_rsSTFT(x, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
end

TFRtot = TFRcurrent ;
        
figure;
imagesc(t,freqs,log1p(abs(TFRcurrent)/5e0)); colormap(1-gray);
axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
drawnow;

if VideoWriting
    myVideo = VideoWriter('../../Results/myVideoFile');
    myVideo.FrameRate = 10 ;
    open(myVideo)
end

k = 1 ;
while n1<Ntot
    t = (n0:n1)/fs ;
    tic;
    x = xtot(n0:n1) ;

    %% Forecasting
    xext = forecasting(x,L,HOP,extK,extM,forecastMethod);
    xx = [x; xext]; % size Nt + L

    %% Update TF Representation
    xxPREC = xx((end-L-round(1.5*Lo*basicTF.hop)):end) ; % Part of the signal to be used to update TFR
    
    switch basicTF.representation
        case 'SST'
            [~, ~, TFRxx, ~] = sqSTFT_C(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 0, 0) ;
        case 'conceFT'
            [~, ~, ~, TFRxx, ~] = ConceFT_sqSTFT_C(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
        case 'STFT'
            [TFRxx, ~] = STFT_C(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10) ;
        case 'RS'
            [~, ~, TFRxx, ~, ~] = ConceFT_rsSTFT(xxPREC, basicTF.fmin, basicTF.fmax, basicTF.df, basicTF.hop, basicTF.win, 1, 10, 1) ;
    end

    TFRcurrent = [TFRcurrent(:,2:(end-Lo)) TFRxx(:,(end-LL-Lo):(end-LL))] ; % sliding sst
    
    dt(k) = toc ;
    
    %% Display
    
    imagesc(t,freqs,log1p(abs(TFRcurrent)/5e0)); %colormap(1-gray); 
    axis xy; % xlabel('Time (s)'); ylabel('Frequency (Hz)');
    drawnow;
    
    n0 = n0 + basicTF.hop ;
    n1 = n0 + Nt - 1 ;
    
    TFRtot = [TFRtot(:,1:(end-Lo)) TFRxx(:,(end-LL-Lo):(end-LL))] ; % whole SST
    k = k+1 ;
    
    if VideoWriting
        frame = getframe(gcf) ;
        writeVideo(myVideo, frame);
    end
    
end

xlabel('Time (s)'); ylabel('Frequency (Hz)');

if VideoWriting
    close(myVideo) ;
end
