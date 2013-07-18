function [ X ] = C(T,N)
%
% C operator performing summation over one dimension of the tensor T \in
% R^(N*N*N)
%

T=reshape(T,N,N,N);
aux=sum(T,1);X=aux(:);
aux=sum(T,2);X=X+aux(:);
aux=sum(T,3);X=X+aux(:);
X=1/3*X;

end

