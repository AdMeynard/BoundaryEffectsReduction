function [Xi,Mu,V] = approxKoopman(X,Y,ordHerm)
% evaluate Koopman modes V, eigen values Mu, and eigen function Psi*Xi
% EDMD (from paper of Williams et al.)

[N,M] = size(X) ;

K = ordHerm^N ; % size of Psi (number of reference functions)
G = zeros(K,K) ;
A = zeros(K,K) ;
for m = 1:M
    psix = HermitN(X(:,m),ordHerm) ;
    psiy = HermitN(Y(:,m),ordHerm) ;
    G = G + (1/M) * (psix' * psix) ;
    A = A + (1/M) * (psix' * psiy) ;
end

KoopMat = pinv(G) * A ;

[Xi,Mu] = eig(KoopMat) ;

B = zeros(K,N) ;
for n = 1:N
    tmp  = ordHerm^(n-1) + 1 ;
    B(tmp,n) = .5 ;
end

V = ( Xi \ B ).' ;
