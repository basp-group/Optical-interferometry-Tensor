function [ T ] = Ct( X,N)
%
%   Ct is the adjoint operator of C
%

X=reshape(X,N,N);
X1=zeros(1,N,N);
X1(1,:,:)=X;


aux=X1(ones(N,1,1),:,:);%equivalent to aux=repmat(X1, [N, 1, 1]);
T=aux;
X1=zeros(N,1,N);
X1(:,1,:)=X;

aux=X1(:,ones(1,N,1),:);%equivalent to aux=repmat(X1, [1, N, 1]);
T=T+aux;
X1=zeros(N,N,1);
X1(:,:,1)=X;

aux=X1(:,:,ones(1,1,N));%equivalent to aux=repmat(X1, [1, 1, N]);
T=T+aux;
T=T/3;
T=T(:);

end

