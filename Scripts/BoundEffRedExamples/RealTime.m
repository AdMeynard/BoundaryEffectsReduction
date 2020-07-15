clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');
addpath('../../Signals/');

%% load THO signal
data = load('PPG.mat');
subT = 2 ; % Sub sampling (pre-processing)
xtot = data.PPG(1:subT:end);
fs = 100/subT ; % sampling frequency
Ntot = length(xtot) ;
ttot = linspace(0,(Ntot-1)/fs,Ntot) ;

VideoWriting = 0 ; % Change to 1 when recording

%% Forecasting parameters
method.name = 'lse' ;
HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round(1.5*L) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

basicTF.hop = round(20/subT) ;
basicTF.win = 2*L+1 ; % window length (in samples)
fmin = 0.5/fs ;
fmax = 4/fs ;
df = 2e-4 ;

freqs = fs * (fmin:df:fmax) ;

Nt = extK+extM+1 ; % segments length
Lo = round(1.5*L/basicTF.hop); % overlap in TF domain
LL = round(L/basicTF.hop); % extension in TF domain
%% Initialization
n0 = 1 ;
n1 = n0+Nt-1 ;

t = ttot(n0:n1) ;
x = xtot(n0:n1) ;
[~, ~, SSTcurrent, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(x, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
SSTtot = SSTcurrent ;
        
figure;
imagesc(t,freqs,log1p(abs(SSTtot)/5e0)); colormap(1-gray);
axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)'); drawnow;

%% Real-time SST
minutes = 2; 
Nmax = minutes*60*fs ;

if VideoWriting
    myVideo = VideoWriter('../../Results/myVideoFile');
    myVideo.FrameRate = 10 ;
    open(myVideo)
end

k = 1 ;
while n1<Nmax
    tic;
    t = ttot(n0:n1) ;
    x = xtot(n0:n1) ;

    %% Forecasting
    xext = forecasting(x,L,HOP,extK,extM,method);
    xxLSE = [x; xext]; % size Nt + L
    
    %% SST
    [~, ~, SSTxx, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE((end-L-round(1.5*Lo*basicTF.hop)):end), fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    
    SSTcurrent = [SSTcurrent(:,2:(end-Lo)) SSTxx(:,(end-LL-Lo):(end-LL))] ; % sliding sst
    
    dt(k) = toc; % time for updating the SST via BoundEffRed

    imagesc(t,freqs,log1p(abs(SSTcurrent)/5e0)); colormap(1-gray); axis xy;
    xlabel('Time (s)'); ylabel('Frequency (Hz)'); drawnow;
    
    n0 = n0 + basicTF.hop ;
    n1 = n0 + Nt - 1 ;
    
    SSTtot = [SSTtot(:,1:(end-Lo)) SSTxx(:,(end-LL-Lo):(end-LL))] ; % whole SST
    k = k+1 ;
    
    if VideoWriting
        frame = getframe(gcf) ;
        writeVideo(myVideo, frame);
    end
    
end

if VideoWriting
    close(myVideo) ;
end

dtmax = basicTF.hop/fs ;
save('../../Results/RealTime','dt','dtmax')
