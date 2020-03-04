function [Xi,mu,phi_end] = approxKoopmanKB(X,Y,sigma2)
% evaluate Koopman modes V, eigen values Mu, and eigen function Psi*Xi
% EDMD (from paper of Hua et al.)

M = size(X,1) ;

Uxy = [X; Y(end,:)] ;

tmp = pdist(Uxy) ;
Uga = exp( -(1/sigma2) * squareform(tmp) ) ;
Uga = Uga(:,1:M);

Ghat = Uga( 1:M, : ) ;
Ahat = Uga( 2:(M+1), : ) ;

[Q, Sigma2] = eig(Ghat) ;
sigmaPINV = 1./sqrt(diag(Sigma2)) ;

Mtmp = bsxfun(@times, Q, sigmaPINV.') ; % <=> Q*pinv(Sigma)
KoopMat = Mtmp.' * Ahat * Mtmp ;
[Vhat, Mu] = eig(KoopMat) ; % eigen values
mu = diag(Mu);

Phixy = Uga * Mtmp * Vhat ;

Xi = pinv(Phixy) * Uxy ;

phi_end = Phixy(end,:);
