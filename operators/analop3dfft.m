function [ y ] = analop3dfft( Tx, n, m, maskB)
%
%
% analop3dfft computes 2D-Fourier Transform along the 3 dimensions of cubic
% tensor T and evaluates it in frequences open in maskB
%
% INPUT
%
%   - Tx: Input tensor
%   - [n, m] dimension of the underlying signal (i.e. Tx \in R^(N*N*N)
%   where N=n*m
%   - maskB mask indicating open frequencies
%


N=n*m;
Tx=reshape(Tx,N,N,N);


y0=analop_fft_6d(Tx,n,m);


maskB=logical(maskB);
 

y=y0(maskB);


end

