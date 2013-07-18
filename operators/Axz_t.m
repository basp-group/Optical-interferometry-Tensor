function [y] = Axz_t( Xhat,xhat,zhat,T )
%
%   Axz_t is the adjoint operator of Axz

n=sqrt(size(xhat,1));m=n;


aux=zeros(size(xhat));
for i=1:size(T,1)
    
    
        aux(T(i,2))=aux(T(i,2))+Xhat(i)*(conj(zhat(T(i,3))*xhat(T(i,1))));
        
end
        
y=n*ifft2(reshape(aux,n,m));       
y=y(:);
end

