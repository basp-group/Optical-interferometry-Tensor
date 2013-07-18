function [x] = Ayz_t( Xhat,yhat,zhat,T )
%
%   Ayz_t is the adjoint operator of Ayz

n=sqrt(size(yhat,1));m=n;


aux=zeros(size(yhat));
for i=1:size(T,1)
    
        aux(T(i,1))=aux(T(i,1))+Xhat(i)*(conj(zhat(T(i,3))*yhat(T(i,2))));
    
    
end
        
x=n*ifft2(reshape(aux,n,m));       
x=x(:);
end

