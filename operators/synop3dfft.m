function [ Tx ] = synop3dfft( y,n,m, maskB)
%
%   synop3dfft is the adjoint operator of analop3dfft
%

N=n*m;

yfull=zeros(N,N,N);
yfull(maskB)=y;

Tx=synop_fft_6d(yfull,n,m);
Tx=Tx(:);
end

