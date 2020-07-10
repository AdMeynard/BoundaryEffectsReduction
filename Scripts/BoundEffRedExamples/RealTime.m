clear all; close all; clc;
addpath(genpath('../../TimeFrequencyScaleRep'));
addpath('../../Algorithm/');
addpath('../../Signals/');

%% load THO signal
data = load('PPG.mat');
xtot = data.PPG;
fs = 100 ; % sampling frequency
Ntot = length(xtot) ;

method.name = 'lse' ;
HOP = 1 ;
extSEC = 5 ; % the extension is of extSEC second
L = round(extSEC*fs) ;
extM = round(1.5*L) ; % dimension of embedding / signals length
extK = round( 2.5*extM ) ;  % number of points to estimate A / size of datasets

basicTF.hop = 10 ;
basicTF.win = 2*L+1 ; % window length (in samples)
fmin = 0/fs ;
fmax = 4/fs ;
df = 2e-5 ;

Nt = extK+extM+1 ; % segments length
Lo = round(1.5*L/basicTF.hop); % overlap in TF domain
LL = round(L/basicTF.hop); % extension in TF domain
%% Initialization
n0 = 1 ;
n1 = n0+Nt-1 ;
x = xtot(n0:n1) ;
[~, ~, SSTtot, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(x, fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;

%n0 = n1-L;
%n1 = n0+Nt-1 ;
        
% figure;
% imagesc(log1p(abs(SSTtot)/5e1)); colormap(1-gray);drawnow;

%% Real-time update
Nmax = 10*60*fs ;

k = 1 ;
while n1<Nmax
    tic;
    x = xtot(n0:n1) ;

    %% Forecasting
    xext = forecasting(x,fs,HOP,extK,extM,extSEC,method);
    xxLSE = [x; xext]; % size Nt + L
%     xxLSE = SigExtension(x,fs,HOP,extK,extM,extSEC,method); 
    
    %% SST
    [~, ~, SSTxx, ~, tfrsqticEXT] = ConceFT_sqSTFT_C(xxLSE((end-L-round(1.5*Lo*basicTF.hop)):end), fmin, fmax,...
            df, basicTF.hop, basicTF.win, 1, 10, 1, 0, 0) ;
    
    SSTtot = [SSTtot(:,1:(end-Lo)) SSTxx(:,(end-LL-Lo):(end-LL))] ;

%     imagesc(log1p(abs(SSTtot)/5e1)); colormap(1-gray);drawnow;
    
    n0 = n0 + basicTF.hop ;
    n1 = n0 + Nt - 1 ;
    
    dt(k) = toc;
    k = k+1 ;
end

dtmax = basicTF.hop/fs ;
save('RealTime','dt','dtmax')
