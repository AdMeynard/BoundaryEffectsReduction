function psix = HermitN(x,K)
% basis function for EDMD (see paper)

N = length(x) ;

psix = hermiteH(0:(K-1), x(1)) ;
for n = 2:N
    tmp = hermiteH(0:(K-1), x(n)) ;
    psix = kron(tmp,psix) ;
end